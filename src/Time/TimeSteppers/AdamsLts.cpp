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

template <typename Op>
LtsCoefficients& add_assign_impl(LtsCoefficients& a, const LtsCoefficients& b,
                                 const Op& op) {
  const auto key_equal = [](const LtsCoefficients::value_type& l,
                            const LtsCoefficients::value_type& r) {
    return get<0>(l) == get<0>(r) and get<1>(l) == get<1>(r);
  };
  const auto key_less = [](const LtsCoefficients::value_type& l,
                           const LtsCoefficients::value_type& r) {
    return get<0>(l) < get<0>(r) or
           (get<0>(l) == get<0>(r) and get<1>(l) < get<1>(r));
  };

  // Two passes: first the common entries, and then the new entries.
  size_t common_entries = 0;
  {
    auto a_it = a.begin();
    auto b_it = b.begin();
    while (a_it != a.end() and b_it != b.end()) {
      if (key_less(*a_it, *b_it)) {
        ++a_it;
      } else if (key_less(*b_it, *a_it)) {
        ++b_it;
      } else {
        ++common_entries;
        get<2>(*a_it) += op(get<2>(*b_it));
        ++a_it;
        ++b_it;
      }
    }
  }
  a.resize(a.size() + b.size() - common_entries);
  {
    auto write = a.rbegin();
    auto a_it = write + static_cast<LtsCoefficients::difference_type>(
                            b.size() - common_entries);
    auto b_it = b.rbegin();
    while (b_it != b.rend()) {
      if (a_it == a.rend() or key_less(*a_it, *b_it)) {
        *write = *b_it;
        get<2>(*write) = op(get<2>(*write));
        ++b_it;
      } else {
        if (key_equal(*a_it, *b_it)) {
          ++b_it;
        }
        *write = *a_it;
        ++a_it;
      }
      ++write;
    }
  }
  return a;
}
}  // namespace

LtsCoefficients& operator+=(LtsCoefficients& a, const LtsCoefficients& b) {
  return add_assign_impl(a, b, [](const double x) { return x; });
}
LtsCoefficients& operator-=(LtsCoefficients& a, const LtsCoefficients& b) {
  return add_assign_impl(a, b, [](const double x) { return -x; });
}

template <typename T>
void apply_coefficients(const gsl::not_null<T*> result,
                        const LtsCoefficients& coefficients,
                        const BoundaryHistoryEvaluator<T>& coupling) {
  for (const auto& term : coefficients) {
    *result += get<2>(term) * *coupling(get<0>(term), get<1>(term));
  }
}

bool operator==(const AdamsScheme& a, const AdamsScheme& b) {
  return a.type == b.type and a.order == b.order;
}
bool operator!=(const AdamsScheme& a, const AdamsScheme& b) {
  return not(a == b);
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

// Collect the ids used for interpolating during a step to `end_time`
// from `times`.
//
// For implicit schemes, the step containing (or ending at) `end_time`
// must have predictor data available.
template <typename TimeType>
OrderVector<TimeStepId> find_relevant_ids2(
    const ConstBoundaryHistoryTimes& times, const TimeType& end_time,
    const AdamsScheme& scheme) {
  OrderVector<TimeStepId> ids{};
  using difference_type = std::iterator_traits<
      ConstBoundaryHistoryTimes::const_iterator>::difference_type;
  const evolution_less<> less{times.front().time_runs_forward()};
  // Can't do a binary search because times are not sorted during self-start.
  auto used_range_end = times.end();
  while (used_range_end != times.begin()) {
    if (less((used_range_end - 1)->step_time(), end_time)) {
      break;
    }
    --used_range_end;
  }
  const auto number_of_past_steps =
      static_cast<difference_type>(scheme.order) -
      (scheme.type == SchemeType::Implicit ? 1 : 0);
  ASSERT(used_range_end - times.begin() >= number_of_past_steps,
         "Insufficient past data.");
  std::copy(used_range_end - number_of_past_steps, used_range_end,
            std::back_inserter(ids));
  ASSERT(scheme.type == SchemeType::Implicit or
             times.number_of_substeps(
                 static_cast<size_t>(used_range_end - times.begin() - 1)) == 1,
         "Taking an explicit step with substep data available.  This is "
         "probably not intended.");
  if (scheme.type == SchemeType::Implicit) {
    const auto last_step =
        static_cast<size_t>(used_range_end - times.begin() - 1);
    ASSERT(times.number_of_substeps(last_step) == 2,
           "Must have substep data for implicit stepping.");
    ids.push_back(times[{last_step, 1}]);
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

// Choose the relevant times from `local` and `remote` for defining
// small steps for `small_step_scheme`, using the specified schemes
// for interpolation on the local and remote sides.
//
// This is the most recent values from the union of the `local` and
// `remote` times with any values that should only be used for
// interpolation removed.
OrderVector<Time> merge_to_small_steps2(const OrderVector<Time>& local,
                                        const OrderVector<Time>& remote,
                                        const evolution_less<Time>& less,
                                        const AdamsScheme& local_scheme,
                                        const AdamsScheme& remote_scheme,
                                        const AdamsScheme& small_step_scheme) {
  ASSERT(small_step_scheme.order <= local.size() and
             small_step_scheme.order <= remote.size(),
         "Insufficient data supplied.");
  OrderVector<Time> small_steps(small_step_scheme.order);
  auto local_it = local.rbegin();
  auto remote_it = remote.rbegin();

  ASSERT(not(local_scheme.type == SchemeType::Implicit and
             remote_scheme.type == SchemeType::Explicit and
             less(*local_it, *remote_it)),
         "Explicit time " << *remote_it << " after implicit " << *local_it);
  ASSERT(not(remote_scheme.type == SchemeType::Implicit and
             local_scheme.type == SchemeType::Explicit and
             less(*remote_it, *local_it)),
         "Explicit time " << *local_it << " after implicit " << *remote_it);

  if (small_step_scheme.type == SchemeType::Explicit) {
    // Don't use implicit interpolation points for an explicit step
    if (local_scheme.type == SchemeType::Implicit) {
      ++local_it;
    }
    if (remote_scheme.type == SchemeType::Implicit) {
      ++remote_it;
    }
  } else {
    if (local_scheme.type == SchemeType::Implicit and
        remote_scheme.type == SchemeType::Implicit) {
      // If both the interpolation schemes are implicit, we will get
      // two times after the small step we are working on, one from
      // each.  One of them (if they are different) belongs to a later
      // small step, and we should ignore it.  If they are the same,
      // ignoring one of them is harmless.
      if (less(*local_it, *remote_it)) {
        ++remote_it;
      } else {
        ++local_it;
      }
    }
  }

  for (auto out = small_steps.rbegin(); out != small_steps.rend(); ++out) {
    if (local_it == local.rend()) {
      ASSERT(remote_it != remote.rend(), "Ran out of data");
      *out = *remote_it;
      ++remote_it;
    } else if (remote_it == remote.rend()) {
      *out = *local_it;
      ++local_it;
    } else {
      *out = std::max(*local_it, *remote_it, less);
      if (*local_it == *out) {
        ++local_it;
      }
      if (*remote_it == *out) {
        ++remote_it;
      }
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

template <typename TimeType>
LtsCoefficients lts_coefficients_for_gts(
    const OrderVector<TimeStepId>& control_ids, const Time& start_time,
    const TimeType& end_time) {
  // The sides are stepping at the same rate, so no LTS is happening
  // at this boundary.
  OrderVector<Time> control_times(control_ids.size());
  std::transform(control_ids.begin(), control_ids.end(), control_times.begin(),
                 exact_substep_time);

  const OrderVector<double> gts_coefficients = adams_coefficients::coefficients(
      control_times.begin(), control_times.end(), start_time, end_time);
  LtsCoefficients lts_coefficients{};
  for (size_t step = 0; step < gts_coefficients.size(); ++step) {
    lts_coefficients.emplace_back(control_ids[step], control_ids[step],
                                  gts_coefficients[step]);
  }
  return lts_coefficients;
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

template <typename TimeType>
LtsCoefficients lts_coefficients2(const ConstBoundaryHistoryTimes& local_times,
                                  const ConstBoundaryHistoryTimes& remote_times,
                                  const Time& start_time,
                                  const TimeType& end_time,
                                  const AdamsScheme& local_scheme,
                                  const AdamsScheme& remote_scheme,
                                  const AdamsScheme& small_step_scheme) {
  if (start_time == end_time) {
    return {};
  }
  const evolution_less<Time> time_less{local_times.front().time_runs_forward()};

  LtsCoefficients step_coefficients{};

  TimeType small_step_end = end_time;
  for (;;) {
    const OrderVector<TimeStepId> local_ids =
        find_relevant_ids2(local_times, small_step_end, local_scheme);
    const OrderVector<TimeStepId> remote_ids =
        find_relevant_ids2(remote_times, small_step_end, remote_scheme);

    // Check is the there is actually local time-stepping happening at
    // this boundary.  Only check for the latest small step, before we
    // have generated any coefficients.
    if (step_coefficients.empty() and small_step_scheme == local_scheme and
        small_step_scheme == remote_scheme and local_ids == remote_ids) {
      return lts_coefficients_for_gts(local_ids, start_time, end_time);
    }

    OrderVector<Time> local_control_times(local_scheme.order);
    std::transform(local_ids.begin(), local_ids.end(),
                   local_control_times.begin(), exact_substep_time);
    OrderVector<Time> remote_control_times(remote_scheme.order);
    std::transform(remote_ids.begin(), remote_ids.end(),
                   remote_control_times.begin(), exact_substep_time);

    const OrderVector<Time> small_step_times = merge_to_small_steps2(
        local_control_times, remote_control_times, time_less, local_scheme,
        remote_scheme, small_step_scheme);
    const Time current_small_step =
        small_step_times[small_step_times.size() -
                         (small_step_scheme.type == SchemeType::Implicit ? 2
                                                                         : 1)];
    ASSERT(not time_less(current_small_step, start_time),
           "Missed step start iterating over small steps.");

    const OrderVector<double> small_step_coefficients =
        adams_coefficients::coefficients(small_step_times.begin(),
                                         small_step_times.end(),
                                         current_small_step, small_step_end);

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
    if (current_small_step == start_time) {
      break;
    }

    if constexpr (std::is_same_v<TimeType, Time>) {
      // We're iterating backwards, so the next step temporally is the
      // one we just did.
      small_step_end = current_small_step;
    } else {
      ERROR(
          "Multiple-small-step dense output is not supported.  Split into an "
          "exact update plus a partial step.");
    }
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
  clean_boundary_history2(local_times, local_times.back().step_time(),
                          integration_order);
  clean_boundary_history2(remote_times, local_times.back().step_time(),
                          integration_order);
}

void clean_boundary_history2(const MutableBoundaryHistoryTimes& times,
                             const Time& first_needed_time,
                             const size_t steps_to_keep) {
  // FIXME take scheme?
  ASSERT(not times.empty(), "Can't clean an empty history.");
  const evolution_less_equal<Time> less_equal{
      times.front().time_runs_forward()};

  // Can't do a binary search because times are not sorted during self-start.
  auto first_needed_step = times.end();
  for (;;) {
    ASSERT(first_needed_step != times.begin(),
           "Trying to clean boundary history for time "
               << first_needed_time << ", but only have data back to "
               << times.front());
    --first_needed_step;
    if (less_equal(first_needed_step->step_time(), first_needed_time)) {
      break;
    }
  }
  const auto number_to_remove = first_needed_step - times.begin() -
                                static_cast<std::ptrdiff_t>(steps_to_keep) + 1;
  ASSERT(number_to_remove >= 0,
         "Insufficient data available during boundary history cleaning.  Have "
             << (first_needed_step - times.begin() + 1)
             << " times available for " << first_needed_time << ", need "
             << steps_to_keep);
  for (int i = 0; i < number_to_remove; ++i) {
    times.pop_front();
  }
  for (size_t i = 0; i < steps_to_keep - 1; ++i) {
    times.clear_substeps(i);
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

#define INSTANTIATE(_, data)                                                 \
  template LtsCoefficients lts_coefficients(                                 \
      const ConstBoundaryHistoryTimes& local_times,                          \
      const ConstBoundaryHistoryTimes& remote_times,                         \
      const TIME_TYPE(data) & end_time, StepType step_type);                 \
  template LtsCoefficients lts_coefficients2(                                \
      const ConstBoundaryHistoryTimes& local_times,                          \
      const ConstBoundaryHistoryTimes& remote_times, const Time& start_time, \
      const TIME_TYPE(data) & end_time, const AdamsScheme& remote_scheme,    \
      const AdamsScheme& local_scheme, const AdamsScheme& small_step_scheme);

GENERATE_INSTANTIATIONS(INSTANTIATE, (Time, ApproximateTime))
#undef INSTANTIATE
#undef TIME_TYPE
}  // namespace TimeSteppers::adams_lts
