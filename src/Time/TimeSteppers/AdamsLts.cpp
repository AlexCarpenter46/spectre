// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/TimeSteppers/AdamsLts.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <type_traits>

#include "DataStructures/MathWrapper.hpp"
#include "NumericalAlgorithms/Interpolation/LagrangePolynomial.hpp"
#include "Time/ApproximateTime.hpp"
#include "Time/BoundaryHistory.hpp"
#include "Time/EvolutionOrdering.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsCoefficients.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace TimeSteppers::adams_lts {
namespace {
template <typename T>
using OrderVector = adams_coefficients::OrderVector<T>;
}  // namespace

template <typename T>
void apply_coefficients(const gsl::not_null<T*> result,
                        const LtsCoefficients& coefficients,
                        const BoundaryHistoryEvaluator<T>& coupling) {
  for (const auto& term : coefficients) {
    *result += get<2>(term) * *coupling(get<0>(term), get<1>(term));
  }
}

namespace {
Time exact_substep_time(const TimeStepId& id) {
  ASSERT(id.substep() == 0 or id.substep() == 1,
         "Implemented Adams-based methods have no more than one substep.");
  if (id.substep() == 0) {
    return id.step_time();
  } else {
    return id.step_time() + id.step_size();
  }
}

// Collect the ids used for interpolating during a step to step number
// `past_end_index` in `times`.
//
// For explicit steps the length of the step doesn't matter, so
// `past_end_index` can be `times.size()`.  For implicit steps, the
// value that would be obtained from `past_end_step` is replaced by
// the predictor from the previous step, so again it can point off the
// end of the vector.
OrderVector<TimeStepId> find_relevant_ids(
    const ConstBoundaryHistoryTimes& times, const size_t past_end_index,
    const size_t number_of_control_times, const bool implicit) {
  OrderVector<TimeStepId> ids{};
  using difference_type = std::iterator_traits<
      ConstBoundaryHistoryTimes::const_iterator>::difference_type;
  const size_t start_offset =
      past_end_index - number_of_control_times + (implicit ? 1 : 0);
  // The usual error triggering this is start_offset wrapping and
  // becoming very large.
  ASSERT(start_offset <= past_end_index, "Insufficient past data.");
  std::copy(times.begin() + static_cast<difference_type>(start_offset),
            times.begin() + static_cast<difference_type>(past_end_index),
            std::back_inserter(ids));
  if (implicit) {
    ASSERT(times.number_of_substeps(past_end_index - 1) == 2,
           "Must have substep data for implicit stepping.");
    ids.push_back(times[{past_end_index - 1, 1}]);
  }
  return ids;
}

// Choose the relevant times from `local` and `remote` for defining
// the small steps, where the Adams coefficients are calculated.  This
// is just the most recent values from the merged list unless
// `drop_partial_implicit` is true, in which case the last value will
// be dropped if it is only to be used for interpolation.
OrderVector<Time> merge_to_small_steps(const OrderVector<Time>& local,
                                       const OrderVector<Time>& remote,
                                       const evolution_less<Time>& less,
                                       const bool drop_partial_implicit) {
  ASSERT(local.size() == remote.size(), "Cannot determine integration order.");
  OrderVector<Time> small_steps(local.size());
  auto local_it = local.rbegin();
  auto remote_it = remote.rbegin();

  if (drop_partial_implicit) {
    // When taking an implicit step, we have data at the step end from
    // both sides, but since the ends of the steps don't necessarily
    // line up we may have even later data on one side.  We don't want
    // to use that as a control point (just for interpolation).
    //
    // For a predictor from an unaligned step, there will, similarly,
    // be a point in the future on the remote side that we want to
    // ignore.
    if (less(*local_it, *remote_it)) {
      // We will be collecting the N largest elements from a list of
      // length N and a list of length N-1, so we have to be sure we
      // won't overrun the shorter source list.  For implicit steps,
      // we work in time intervals between the passed times, so the
      // last interval on each side must overlap or we've been passed
      // useless information.  This guarantees that the first time
      // selected will be taken from the longer list only, so we won't
      // overrun the shorter list.  This is ASSERTed below.
      //
      // That first selection still examines the largest element on
      // the shorter list even though it will not be chosen, however,
      // so we must handle the case where N=1 specially.  This can
      // only happen for a first-order predictor for an unaligned
      // step.
      if (small_steps.size() == 1) {
        return local;
      }
      ++remote_it;
      ASSERT(less(*remote_it, *local_it), "Received extra remote step");
    } else {
      ASSERT(small_steps.size() != 1,
             "First order steps can only be implicit on the remote side.");
      ++local_it;
      ASSERT(less(*local_it, *remote_it), "Received extra local step");
    }
  }
  for (auto out = small_steps.rbegin(); out != small_steps.rend(); ++out) {
    *out = std::max(*local_it, *remote_it, less);
    if (*local_it == *out) {
      ++local_it;
    }
    if (*remote_it == *out) {
      ++remote_it;
    }
  }
  return small_steps;
}

// Evaluate the Lagrange interpolating polynomials with the given
// `control_times` at `time`.  The returned vector contains the values
// of all the Lagrange polynomials.
template <typename TimeType>
OrderVector<double> interpolation_coefficients(
    const OrderVector<Time>& control_times, const TimeType& time) {
  if constexpr (std::is_same_v<TimeType, Time>) {
    // Skip the Lagrange polynomial calculations if we are evaluating
    // at a control time.  This should be common.
    for (size_t i = 0; i < control_times.size(); ++i) {
      if (control_times[i] == time) {
        OrderVector<double> coefficients(control_times.size(), 0.0);
        coefficients[i] = 1.0;
        return coefficients;
      }
    }
  }

  OrderVector<double> control_times_fp(control_times.size());
  std::transform(control_times.begin(), control_times.end(),
                 control_times_fp.begin(),
                 [](const Time& t) { return t.value(); });

  OrderVector<double> coefficients{};
  for (size_t i = 0; i < control_times.size(); ++i) {
    coefficients.push_back(lagrange_polynomial(
        i, time.value(), control_times_fp.begin(), control_times_fp.end()));
  }
  return coefficients;
}
}  // namespace

template <typename TimeType>
LtsCoefficients lts_coefficients(const ConstBoundaryHistoryTimes& local_times,
                                 const ConstBoundaryHistoryTimes& remote_times,
                                 const TimeType& end_time,
                                 const StepType step_type) {
  ASSERT((std::is_same_v<TimeType, Time>) or step_type != StepType::Predictor,
         "Can't do dense output on predictor stage.");
  const size_t order = local_times.integration_order(local_times.size() - 1) -
                       (step_type == StepType::Predictor ? 1 : 0);
  ASSERT(order > 0, "Zeroth-order LTS not supported.");
  ASSERT(step_type != StepType::Corrector or order > 1,
         "First order predictor-corrector not supported.");
  const evolution_less<Time> time_less{local_times.front().time_runs_forward()};

  const bool implicit = step_type == StepType::Corrector;
  // The predictor uses implicit data on the remote side if the sides
  // aren't aligned at the step start, because the remote predictor
  // evaluation is sequenced before the interior of its step.
  const bool remote_implicit =
      implicit or (step_type == StepType::Predictor and
                   local_times.back() != remote_times.back());

  const Time step_start = local_times.back().step_time();

  const OrderVector<TimeStepId> local_ids =
      find_relevant_ids(local_times, local_times.size(), order, implicit);
  OrderVector<Time> local_control_times(order);
  std::transform(local_ids.begin(), local_ids.end(),
                 local_control_times.begin(), exact_substep_time);

  // These are updated in the loop below.
  size_t remote_index = remote_times.size();
  OrderVector<TimeStepId> remote_ids =
      find_relevant_ids(remote_times, remote_index, order, remote_implicit);

  if (local_ids == remote_ids) {
    // The sides are stepping at the same rate, so no LTS is happening
    // at this boundary.
    const OrderVector<double> gts_coefficients =
        adams_coefficients::coefficients(local_control_times.begin(),
                                         local_control_times.end(), step_start,
                                         end_time);
    LtsCoefficients lts_coefficients{};
    for (size_t step = 0; step < order; ++step) {
      lts_coefficients.emplace_back(local_ids[step], local_ids[step],
                                    gts_coefficients[step]);
    }
    return lts_coefficients;
  }

  LtsCoefficients step_coefficients{};

  Time next_small_step{};
  for (;;) {
    OrderVector<Time> remote_control_times(order);
    std::transform(remote_ids.begin(), remote_ids.end(),
                   remote_control_times.begin(), exact_substep_time);

    // implicit implies remote_implicit, so we don't need to check it.
    const OrderVector<Time> small_step_times = merge_to_small_steps(
        local_control_times, remote_control_times, time_less, remote_implicit);
    const Time current_small_step =
        small_step_times[small_step_times.size() - (implicit ? 2 : 1)];
    const OrderVector<double> small_step_coefficients =
        step_coefficients.empty()
            ? adams_coefficients::coefficients(small_step_times.begin(),
                                               small_step_times.end(),
                                               current_small_step, end_time)
            : adams_coefficients::coefficients(
                  small_step_times.begin(), small_step_times.end(),
                  current_small_step, next_small_step);

    for (size_t contributing_small_step = 0;
         contributing_small_step < small_step_times.size();
         ++contributing_small_step) {
      const OrderVector<double> local_interpolation_coefficients =
          interpolation_coefficients(local_control_times,
                                     small_step_times[contributing_small_step]);
      const OrderVector<double> remote_interpolation_coefficients =
          interpolation_coefficients(remote_control_times,
                                     small_step_times[contributing_small_step]);
      for (size_t local_step_index = 0;
           local_step_index < local_interpolation_coefficients.size();
           ++local_step_index) {
        if (local_interpolation_coefficients[local_step_index] == 0.0) {
          continue;
        }
        for (size_t remote_step_index = 0;
             remote_step_index < remote_interpolation_coefficients.size();
             ++remote_step_index) {
          if (remote_interpolation_coefficients[remote_step_index] == 0.0) {
            continue;
          }
          step_coefficients.emplace_back(
              local_ids[local_step_index], remote_ids[remote_step_index],
              small_step_coefficients[contributing_small_step] *
                  local_interpolation_coefficients[local_step_index] *
                  remote_interpolation_coefficients[remote_step_index]);
        }
      }
    }
    if (current_small_step == step_start) {
      break;
    }
    ASSERT(not time_less(current_small_step, step_start),
           "Missed step start iterating over small steps.");
    ASSERT(step_type != StepType::Predictor,
           "Predictor should only have one small step.");

    --remote_index;
    remote_ids =
        find_relevant_ids(remote_times, remote_index, order, remote_implicit);
    // We're iterating backwards, so the next step temporally is the
    // one we just did.
    next_small_step = current_small_step;
  }

  // Combine duplicate entries
  ASSERT(not step_coefficients.empty(), "Generated no coefficients");
  alg::sort(step_coefficients);
  auto unique_entry = step_coefficients.begin();
  for (auto generated_entry = std::next(step_coefficients.begin());
       generated_entry != step_coefficients.end();
       ++generated_entry) {
    if (get<0>(*generated_entry) == get<0>(*unique_entry) and
        get<1>(*generated_entry) == get<1>(*unique_entry)) {
      get<2>(*unique_entry) += get<2>(*generated_entry);
    } else {
      *++unique_entry = *generated_entry;
    }
  }
  step_coefficients.erase(std::next(unique_entry), step_coefficients.end());
  return step_coefficients;
}

void clean_boundary_history(const MutableBoundaryHistoryTimes& local_times,
                            const MutableBoundaryHistoryTimes& remote_times,
                            const size_t integration_order) {
  ASSERT(local_times.number_of_substeps(local_times.size() - 1) == 1,
         "Last local step should not have substep data during cleaning.");

  const auto signed_order =
      static_cast<typename decltype(local_times.end())::difference_type>(
          integration_order);

  ASSERT(local_times.size() >= integration_order,
         "Insufficient data to take an order-" << integration_order
         << " step.  Have " << local_times.size() << " times, need "
         << integration_order);
  while (local_times.size() > integration_order) {
    local_times.pop_front();
  }
  for (size_t step = 0; step < integration_order - 1; ++step) {
    local_times.clear_substeps(step);
  }

  ASSERT(remote_times.size() >= integration_order,
         "Insufficient data to take an order-" << integration_order
         << " step.  Have " << remote_times.size() << " times, need "
         << integration_order);
  if (std::equal(local_times.begin(), local_times.end(),
                 remote_times.end() - signed_order)) {
    // GTS
    while (remote_times.size() > integration_order) {
      remote_times.pop_front();
    }
    for (size_t step = 0; step < integration_order; ++step) {
      remote_times.clear_substeps(step);
    }
  } else {
    const auto remote_step_after_step_start = std::upper_bound(
        remote_times.begin(), remote_times.end(), local_times.back());
    ASSERT(remote_step_after_step_start - remote_times.begin() >= signed_order,
           "Insufficient data to take an order-" << integration_order
           << " step.  Have "
           << remote_step_after_step_start - remote_times.begin()
           << " times before the step, need " << integration_order);
    // The pop_front() calls invalidate remote_step_after_step_start.
    const TimeStepId first_needed =
        *(remote_step_after_step_start - signed_order);
    while (remote_times.front() != first_needed) {
      remote_times.pop_front();
    }
    const evolution_less_equal<Time> less_equal{
        local_times.front().time_runs_forward()};
    for (size_t step = 0; step < integration_order - 1; ++step) {
      if (remote_times.number_of_substeps(step) == 2 and
          less_equal(exact_substep_time(remote_times[{step, 1}]),
                     local_times.back().step_time())) {
        remote_times.clear_substeps(step);
      }
    }
  }
}

#define MATH_WRAPPER_TYPE(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                          \
  template void apply_coefficients(                   \
      gsl::not_null<MATH_WRAPPER_TYPE(data)*> result, \
      const LtsCoefficients& coefficients,            \
      const BoundaryHistoryEvaluator<MATH_WRAPPER_TYPE(data)>& coupling);

GENERATE_INSTANTIATIONS(INSTANTIATE, (MATH_WRAPPER_TYPES))
#undef INSTANTIATE
#undef MATH_WRAPPER_TYPE

#define TIME_TYPE(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                         \
  template LtsCoefficients lts_coefficients(         \
      const ConstBoundaryHistoryTimes& local_times,  \
      const ConstBoundaryHistoryTimes& remote_times, \
      const TIME_TYPE(data) & end_time, StepType step_type);

GENERATE_INSTANTIATIONS(INSTANTIATE, (Time, ApproximateTime))
#undef INSTANTIATE
#undef TIME_TYPE
}  // namespace TimeSteppers::adams_lts
