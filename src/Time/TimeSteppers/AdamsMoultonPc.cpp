// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/TimeSteppers/AdamsMoultonPc.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <pup.h>

#include "Time/ApproximateTime.hpp"
#include "Time/EvolutionOrdering.hpp"
#include "Time/History.hpp"
#include "Time/SelfStart.hpp"
#include "Time/TimeSteppers/AdamsCoefficients.hpp"
#include "Time/TimeSteppers/AdamsLts.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

// FIXME make both available
#define MONOTONIC

namespace TimeSteppers {

// Don't include AdamsCoefficients.hpp in the header just to get one
// constant.
static_assert(adams_coefficients::maximum_order ==
              AdamsMoultonPc::maximum_order);

namespace {
// FIXME copied from AdamsLts.cpp
Time exact_substep_time(const TimeStepId& id) {
  ASSERT(id.substep() == 0 or id.substep() == 1,
         "Implemented Adams-based methods have no more than one substep.");
  if (id.substep() == 0) {
    return id.step_time();
  } else {
    return id.step_time() + id.step_size();
  }
}

template <typename T>
void clean_history(const MutableUntypedHistory<T>& history) {
#ifndef MONOTONIC
  ASSERT(history.integration_order() > 1, "Cannot run below second order.");
  const auto required_points = history.integration_order() - 1;
  ASSERT(history.size() >= required_points,
         "Insufficient data to take an order-" << history.integration_order()
         << " step.  Have " << history.size() << " times, need "
         << required_points);
  if (history.at_step_start()) {
    history.clear_substeps();
    while (history.size() > required_points) {
      history.pop_front();
    }
    if (history.size() > 1) {
      history.discard_value(history[history.size() - 2].time_step_id);
    }
  } else {
    history.discard_value(history.substeps().back().time_step_id);
  }
#endif
}

template <typename T>
void clean_history2(const MutableUntypedHistory<T>& history) {
#ifdef MONOTONIC
  ASSERT(history.integration_order() > 1, "Cannot run below second order.");
  // FIXME Steady-state only needs -2, but self-start needs -1.  Make
  // conditional?
  //const auto required_points = history.integration_order() - 2;
  const auto required_points = history.integration_order() - 1;
  ASSERT(history.size() >= required_points,
         "Insufficient data to take an order-" << history.integration_order()
         << " step.  Have " << history.size() << " times, need "
         << required_points);
  if (not history.at_step_start()) {
    history.clear_substeps();
    while (history.size() > required_points) {
      history.pop_front();
    }
    if (not history.empty()) {
      history.discard_value(history.back().time_step_id);
    }
  }
#endif
}

template <typename T, typename TimeType>
void update_u_common(const gsl::not_null<T*> u,
                     const ConstUntypedHistory<T>& history,
                     const TimeType& step_end, const size_t method_order,
                     const bool corrector) {
  ASSERT(history.size() >= method_order - 1, "Insufficient history");
  // Pass in whether to run the predictor or corrector even though we
  // can compute it as a sanity check.
  ASSERT(corrector != history.substeps().empty(),
         "Applying predictor or corrector when expecting the other.");
  ASSERT(corrector != history.at_step_start(), "Unexpected new data");

  const auto used_history_begin =
      history.end() -
      static_cast<typename ConstUntypedHistory<T>::difference_type>(
          method_order - 1);
  adams_coefficients::OrderVector<Time> control_times{};
  std::transform(used_history_begin, history.end(),
                 std::back_inserter(control_times),
                 [](const auto& r) { return r.time_step_id.step_time(); });
  if (corrector) {
    control_times.push_back(
        history.back().time_step_id.step_time() +
        history.substeps().front().time_step_id.step_size());
  }
  const auto coefficients = adams_coefficients::coefficients(
      control_times.begin(), control_times.end(),
      history.back().time_step_id.step_time(), step_end);

  auto coefficient = coefficients.begin();
  for (auto history_entry = used_history_begin;
       history_entry != history.end();
       ++history_entry, ++coefficient) {
    *u += *coefficient * history_entry->derivative;
  }
  if (corrector) {
    *u += coefficients.back() * history.substeps().front().derivative;
  }
}
}  // namespace

AdamsMoultonPc::AdamsMoultonPc(const size_t order) : order_(order) {
  ASSERT(order >= minimum_order and order <= maximum_order,
         "Invalid order: " << order);
}

size_t AdamsMoultonPc::order() const { return order_; }

size_t AdamsMoultonPc::error_estimate_order() const { return order_ - 1; }

uint64_t AdamsMoultonPc::number_of_substeps() const { return 2; }

uint64_t AdamsMoultonPc::number_of_substeps_for_error() const {
  return number_of_substeps();
}

size_t AdamsMoultonPc::number_of_past_steps() const { return order_ - 2; }

double AdamsMoultonPc::stable_step() const {
  switch (order_) {
    case 2:
      return 1.0;
    case 3:
      return 0.981297;
    case 4:
      return 0.794227;
    case 5:
      return 0.612340;
    case 6:
      return 0.464542;
    case 7:
      return 0.350596;
    case 8:
      return 0.264373;
    default:
      ERROR("Bad order");
  }
}

TimeStepId AdamsMoultonPc::next_time_id(const TimeStepId& current_id,
                                        const TimeDelta& time_step) const {
  switch (current_id.substep()) {
    case 0:
      return current_id.next_substep(time_step, 1.0);
    case 1:
      return current_id.next_step(time_step);
    default:
      ERROR("Bad id: " << current_id);
  }
}

TimeStepId AdamsMoultonPc::next_time_id_for_error(
    const TimeStepId& current_id, const TimeDelta& time_step) const {
  return next_time_id(current_id, time_step);
}

bool AdamsMoultonPc::neighbor_data_required(
    const TimeStepId& next_substep_id,
    const TimeStepId& neighbor_data_id) const {
  const evolution_less<Time> before{neighbor_data_id.time_runs_forward()};

  // FIXME self-start - can this be combined somehow?
  if (neighbor_data_id.slab_number() != next_substep_id.slab_number()) {
    return neighbor_data_id.slab_number() < next_substep_id.slab_number();
  }

#ifdef MONOTONIC
  const auto next_time = exact_substep_time(next_substep_id);
  const auto neighbor_time = exact_substep_time(neighbor_data_id);
  return before(neighbor_time, next_time) or
         (neighbor_time == next_time and neighbor_data_id.substep() == 1 and
          next_substep_id.substep() == 0);
#else
  if (next_substep_id.substep() == 1) {
    // predictor
    return before(neighbor_data_id.step_time(), next_substep_id.step_time()) or
           (neighbor_data_id.step_time() == next_substep_id.step_time() and
            neighbor_data_id.substep() == 0);
  } else {
    // corrector
    return before(neighbor_data_id.step_time(), next_substep_id.step_time());
  }
#endif
}

bool AdamsMoultonPc::neighbor_data_required(
    const double dense_output_time, const TimeStepId& neighbor_data_id) const {
  const evolution_less<double> before{neighbor_data_id.time_runs_forward()};
#ifdef MONOTONIC
  return before(exact_substep_time(neighbor_data_id).value(),
                dense_output_time) or
         (neighbor_data_id.substep() == 1 and
          exact_substep_time(neighbor_data_id).value() == dense_output_time);
#else
  return before(neighbor_data_id.step_time().value(), dense_output_time);
#endif
}

void AdamsMoultonPc::pup(PUP::er& p) {
  TimeStepper::pup(p);
  p | order_;
}

template <typename T>
void AdamsMoultonPc::update_u_impl(const gsl::not_null<T*> u,
                                   const MutableUntypedHistory<T>& history,
                                   const TimeDelta& time_step) const {
  clean_history(history);
  const Time next_time = history.back().time_step_id.step_time() + time_step;
  *u = *history.back().value;
  update_u_common(u, history, next_time, history.integration_order(),
                  not history.at_step_start());
  // FIXME can we clean at end?  document step redoing guarantees
  clean_history2(history);
}

template <typename T>
bool AdamsMoultonPc::update_u_impl(const gsl::not_null<T*> u,
                                   const gsl::not_null<T*> u_error,
                                   const MutableUntypedHistory<T>& history,
                                   const TimeDelta& time_step) const {
  clean_history(history);
  const bool predictor = history.at_step_start();
  const Time next_time = history.back().time_step_id.step_time() + time_step;
  *u = *history.back().value;
  update_u_common(u, history, next_time, history.integration_order(),
                  not predictor);
  if (predictor) {
    return false;
  }
  *u_error = *history.back().value;
  update_u_common(u_error, history, next_time, history.integration_order() - 1,
                  true);
  *u_error = *u - *u_error;
  clean_history2(history);
  return true;
}

template <typename T>
bool AdamsMoultonPc::dense_update_u_impl(const gsl::not_null<T*> u,
                                         const ConstUntypedHistory<T>& history,
                                         const double time) const {
  // Special case required to handle the initial time.
  if (time == history.back().time_step_id.step_time().value()) {
    return true;
  }
#ifdef MONOTONIC
  if (not history.at_step_start()) {
    return false;
  }
  update_u_common(u, history, ApproximateTime{time},
                  history.integration_order(), false);
#else
  if (history.at_step_start()) {
    return false;
  }
  update_u_common(u, history, ApproximateTime{time},
                  history.integration_order(), true);
#endif
  return true;
}

template <typename T>
bool AdamsMoultonPc::can_change_step_size_impl(
    const TimeStepId& time_id, const ConstUntypedHistory<T>& history) const {
  // We need to prevent the next step from occurring at the same time
  // as one already in the history.  The self-start code ensures this
  // can't happen during self-start, and it clearly can't happen
  // during normal evolution where the steps are monotonic, but during
  // the transition between them we have to worry about a step being
  // placed on a self-start time.  The self-start algorithm guarantees
  // the final state is safe for constant-time-step evolution, so we
  // just force that until we've passed all the self-start times.
  const evolution_less_equal<Time> less_equal{time_id.time_runs_forward()};
  return not ::SelfStart::is_self_starting(time_id) and
         time_id.substep() == 0 and
         alg::all_of(history, [&](const auto& record) {
           return less_equal(record.time_step_id.step_time(),
                             time_id.step_time());
         });
}

template <typename T>
void AdamsMoultonPc::add_boundary_delta_impl(
    const gsl::not_null<T*> result,
    const TimeSteppers::MutableBoundaryHistoryTimes& local_times,
    const TimeSteppers::MutableBoundaryHistoryTimes& remote_times,
    const TimeSteppers::BoundaryHistoryEvaluator<T>& coupling,
    const TimeDelta& time_step) const {
  ASSERT(not local_times.empty(), "No local data provided.");
  ASSERT(not remote_times.empty(), "No remote data provided.");
  const auto current_order =
      local_times.integration_order(local_times.size() - 1);
#ifdef MONOTONIC
  const adams_lts::AdamsScheme predictor_scheme{adams_lts::SchemeType::Explicit,
                                                current_order - 1};
  const adams_lts::AdamsScheme corrector_scheme{adams_lts::SchemeType::Implicit,
                                                current_order};
  const auto step_start = local_times.back().step_time();
  const auto step_end = step_start + time_step;
  const auto small_step_start =
      std::max(local_times.back(), remote_times.back()).step_time();
  const auto synchronization_time =
      std::min(local_times.back(), remote_times.back()).step_time();
  const auto is_synchronization_time = [&](const TimeStepId& id) {
    return id.step_time() == synchronization_time;
  };
  ASSERT(alg::any_of(local_times, is_synchronization_time) and
             alg::any_of(remote_times, is_synchronization_time),
         "Only nested step patterns (N:1) are supported.");

  if (local_times.number_of_substeps(local_times.size() - 1) == 1) {
    // Predictor
    adams_lts::clean_boundary_history2(local_times, synchronization_time,
                                       current_order - 1);
    adams_lts::clean_boundary_history2(remote_times, synchronization_time,
                                       current_order - 1);

    auto lts_coefficients = adams_lts::lts_coefficients2(
        local_times, remote_times, small_step_start, step_end, predictor_scheme,
        predictor_scheme, predictor_scheme);
    if (not is_synchronization_time(remote_times.back())) {
      lts_coefficients += adams_lts::lts_coefficients2(
          local_times, remote_times, synchronization_time, small_step_start,
          predictor_scheme, corrector_scheme, corrector_scheme);
    }
    adams_lts::apply_coefficients(result, lts_coefficients, coupling);
  } else {
    // Corrector
    if (remote_times.number_of_substeps(remote_times.size() - 1) == 2) {
      // Aligned corrector
      ASSERT(exact_substep_time(remote_times[{remote_times.size() - 1, 1}]) ==
             step_end,
             "Have remote substep data, but it isn't aligned with the local "
             "data.");
      const auto lts_coefficients =
          adams_lts::lts_coefficients2(local_times, remote_times,
                                       synchronization_time, step_end,
                                       corrector_scheme, corrector_scheme,
                                       corrector_scheme) -
          adams_lts::lts_coefficients2(local_times, remote_times,
                                       synchronization_time, step_start,
                                       corrector_scheme, predictor_scheme,
                                       corrector_scheme);
      adams_lts::apply_coefficients(result, lts_coefficients, coupling);
    } else {
      // Unaligned corrector
      ASSERT(step_start == small_step_start,
             "Trying to take unaligned step, but remote side is smaller.");
      const auto lts_coefficients = adams_lts::lts_coefficients2(
          local_times, remote_times, step_start, step_end, corrector_scheme,
          predictor_scheme, corrector_scheme);
      adams_lts::apply_coefficients(result, lts_coefficients, coupling);
    }
  }
#else
  adams_lts::AdamsScheme scheme{adams_lts::SchemeType::Implicit, current_order};
  auto remote_scheme = scheme;

  if (local_times.number_of_substeps(local_times.size() - 1) == 1) {
    // Predictor
    adams_lts::clean_boundary_history(local_times, remote_times,
                                      current_order - 1);

    scheme = {adams_lts::SchemeType::Explicit, current_order - 1};
    ASSERT(remote_times.back() <= local_times.back(),
           "Unexpected remote values available.");
    // If the sides are not aligned, we use the predictor data
    // available from the neighbor.  If they are, that data has not
    // been received.
    if (remote_times.back() == local_times.back()) {
      remote_scheme = scheme;
    } else {
      // FIXME do we want to do this?  If so, clearer ways to write code.
      // Change the interpolation order to the predictor order.
      remote_scheme.order = current_order - 1;
    }
  }

  const auto lts_coefficients = adams_lts::lts_coefficients2(
      local_times, remote_times, local_times.back().step_time(),
      local_times.back().step_time() + time_step, scheme, remote_scheme,
      scheme);
  adams_lts::apply_coefficients(result, lts_coefficients, coupling);
#endif
}

template <typename T>
void AdamsMoultonPc::boundary_dense_output_impl(
    const gsl::not_null<T*> result,
    const TimeSteppers::ConstBoundaryHistoryTimes& local_times,
    const TimeSteppers::ConstBoundaryHistoryTimes& remote_times,
    const TimeSteppers::BoundaryHistoryEvaluator<T>& coupling,
    const double time) const {
#ifdef MONOTONIC
  ASSERT(local_times.number_of_substeps(local_times.size() - 1) == 1,
         "Dense output must be done before predictor evaluation.");

  const auto current_order =
      local_times.integration_order(local_times.size() - 1);
  const adams_lts::AdamsScheme predictor_scheme{adams_lts::SchemeType::Explicit,
                                                current_order - 1};
  const adams_lts::AdamsScheme corrector_scheme{adams_lts::SchemeType::Implicit,
                                                current_order};
  const auto small_step_start =
      std::max(local_times.back(), remote_times.back()).step_time();
  const auto synchronization_time =
      std::min(local_times.back(), remote_times.back()).step_time();
  const auto is_synchronization_time = [&](const TimeStepId& id) {
    return id.step_time() == synchronization_time;
  };
  ASSERT(alg::any_of(local_times, is_synchronization_time) and
             alg::any_of(remote_times, is_synchronization_time),
         "Only nested step patterns (N:1) are supported.");

  auto lts_coefficients = adams_lts::lts_coefficients2(
      local_times, remote_times, small_step_start, ApproximateTime{time},
      predictor_scheme, predictor_scheme, predictor_scheme);
  if (not is_synchronization_time(remote_times.back())) {
    lts_coefficients += adams_lts::lts_coefficients2(
        local_times, remote_times, synchronization_time, small_step_start,
        predictor_scheme, corrector_scheme, corrector_scheme);
  }
  adams_lts::apply_coefficients(result, lts_coefficients, coupling);
#else
  if (local_times.back().step_time().value() == time) {
    // Nothing to do.  The requested time is the start of the step,
    // which is the input value of `result`.
    return;
  }
  ASSERT(local_times.number_of_substeps(local_times.size() - 1) == 2,
         "Dense output must be done after predictor evaluation.");

  const auto current_order =
      local_times.integration_order(local_times.size() - 1);
  const adams_lts::AdamsScheme scheme{adams_lts::SchemeType::Implicit,
                                      current_order};
  const auto small_step_start =
      std::max(local_times.back(), remote_times.back()).step_time();
  const auto lts_coefficients =
      adams_lts::lts_coefficients2(local_times, remote_times,
                                   local_times.back().step_time(),
                                   small_step_start, scheme, scheme, scheme) +
      adams_lts::lts_coefficients2(local_times, remote_times, small_step_start,
                                   ApproximateTime{time}, scheme, scheme,
                                   scheme);
  adams_lts::apply_coefficients(result, lts_coefficients, coupling);
#endif
}

bool operator==(const AdamsMoultonPc& lhs, const AdamsMoultonPc& rhs) {
  return lhs.order() == rhs.order();
}

bool operator!=(const AdamsMoultonPc& lhs, const AdamsMoultonPc& rhs) {
  return not(lhs == rhs);
}

TIME_STEPPER_DEFINE_OVERLOADS(AdamsMoultonPc)
LTS_TIME_STEPPER_DEFINE_OVERLOADS(AdamsMoultonPc)
}  // namespace TimeSteppers

PUP::able::PUP_ID TimeSteppers::AdamsMoultonPc::my_PUP_ID = 0;  // NOLINT
