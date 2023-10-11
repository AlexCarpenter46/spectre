// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Time/TimeSteppers/LtsHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>

#include "Helpers/Time/TimeSteppers/TimeStepperTestUtils.hpp"
#include "Time/BoundaryHistory.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace TimeStepperTestUtils::lts {
namespace {
// Initialize a BoundaryHistory for each side of an interface with
// constant steps if the integrator requires such data.  No entry is
// added for the initial time.
template <
    typename CouplingResult, typename Func1, typename Func2,
    typename Data1 = std::decay_t<std::invoke_result_t<const Func1&, double>>,
    typename Data2 = std::decay_t<std::invoke_result_t<const Func2&, double>>>
auto initialize_history(const size_t number_of_past_values, const size_t order,
                        const double initial_time, const double step_size,
                        const Func1& data_func1, const Func2& data_func2)
    -> std::pair<TimeSteppers::BoundaryHistory<Data1, Data2, CouplingResult>,
                 TimeSteppers::BoundaryHistory<Data2, Data1, CouplingResult>> {
  std::pair<TimeSteppers::BoundaryHistory<Data1, Data2, CouplingResult>,
            TimeSteppers::BoundaryHistory<Data2, Data1, CouplingResult>>
      result;

  const bool time_runs_forward = step_size > 0.0;
  // We use one step per slab.
  auto slab = time_runs_forward
                  ? Slab::with_duration_to_end(initial_time, step_size)
                  : Slab::with_duration_from_start(initial_time, -step_size);
  const auto first_slab_number = -static_cast<int64_t>(number_of_past_values);
  for (int64_t slab_num = -1; slab_num >= first_slab_number; --slab_num) {
    const TimeStepId id(time_runs_forward, slab_num,
                        time_runs_forward ? slab.start() : slab.end());
    const double time = id.step_time().value();
    const auto data1 = data_func1(time);
    const auto data2 = data_func2(time);
    result.first.local().insert_initial(id, order, data1);
    result.second.remote().insert_initial(id, order, data1);
    result.second.local().insert_initial(id, order, data2);
    result.first.remote().insert_initial(id, order, data2);
    slab = time_runs_forward ? slab.retreat() : slab.advance();
  }
  return result;
}

template <typename Rhs>
struct SideData {
  TimeStepId id;
  TimeStepId next_id;
  double step_value;
  double substep_value;
  TimeSteppers::BoundaryHistory<double, double, double> history;
  std::deque<std::pair<TimeStepId, double>> messages;
  Rhs rhs;
  Rational (*pattern)(const Rational&);
};

// FIXME dense output
template <typename Pattern>
double lts_error(const LtsTimeStepper& stepper, const bool time_runs_forward,
                 const int32_t pattern_repeats,
                 const std::optional<double> dense_output_fraction) {
  ASSERT(not dense_output_fraction.has_value() or
             (*dense_output_fraction >= 0.0 and *dense_output_fraction <= 1.0),
         "Bad dense output fraction: " << dense_output_fraction);

  // Test system:
  // dx/dt = x y
  // dy/dt = - x y
  //
  // Solution:
  // x = c / [1 + exp(-c (t - d))]
  // y = c - x
  const auto rhs_x = [](const double x, const double y) { return x * y; };
  const auto rhs_y = [](const double x, const double y) { return -x * y; };

  // Arbitrary
  double conserved_sum = 0.7;
  double offset = 0.4;
  double t_init = 0.6;
  double t_final = 6.3;
  if (not time_runs_forward) {
    conserved_sum *= -1.0;
    offset *= -1.0;
    t_init *= -1.0;
    t_final *= -1.0;
  }

  const auto expected_x = [&conserved_sum, &offset](const double t) {
    return conserved_sum / (1.0 + exp(-conserved_sum * (t - offset)));
  };
  const auto expected_y = [&conserved_sum, &expected_x](const double t) {
    return conserved_sum - expected_x(t);
  };

  const auto slab =
      time_runs_forward ? Slab(t_init, t_final) : Slab(t_final, t_init);
  const TimeDelta pattern_size =
      (time_runs_forward ? 1 : -1) * slab.duration() / pattern_repeats;
  const TimeStepId initial_id(time_runs_forward, 0,
                              time_runs_forward ? slab.start() : slab.end());
  const TimeStepId final_id(time_runs_forward, 0,
                            time_runs_forward ? slab.end() : slab.start());

  auto histories = initialize_history<double>(
      stepper.number_of_past_steps(), stepper.order(), t_init,
      (t_final - t_init) / pattern_repeats, expected_x, expected_y);

  histories.first.local().insert(initial_id, stepper.order(),
                                 expected_x(t_init));
  histories.second.remote().insert(initial_id, stepper.order(),
                                   expected_x(t_init));
  histories.second.local().insert(initial_id, stepper.order(),
                                  expected_y(t_init));
  histories.first.remote().insert(initial_id, stepper.order(),
                                  expected_y(t_init));

  SideData<decltype(&rhs_x)> data_x{
      initial_id,
      stepper.next_time_id(initial_id, Pattern::x(0) * pattern_size),
      expected_x(t_init),
      expected_x(t_init),
      std::move(histories.first),
      {},
      &rhs_x,
      &Pattern::x};

  SideData<decltype(&rhs_y)> data_y{
      initial_id,
      stepper.next_time_id(initial_id, Pattern::y(0) * pattern_size),
      expected_y(t_init),
      expected_y(t_init),
      std::move(histories.second),
      {},
      &rhs_y,
      &Pattern::y};

  // FIXME make function?
  const auto take_steps = [&final_id, &pattern_size, &stepper](
                              const auto local_data, const auto remote_data) {
    while (not stepper.neighbor_data_required(local_data->next_id,
                                              remote_data->next_id)) {
      auto& local_messages = local_data->messages;
      while (not local_messages.empty() and
             stepper.neighbor_data_required(local_data->next_id,
                                            local_messages.front().first)) {
        local_data->history.remote().insert(local_messages.front().first,
                                            stepper.order(),
                                            local_messages.front().second);
        local_messages.pop_front();
      }

      const auto pattern_fraction_after =
          [&pattern_size](const TimeStepId& id) {
            Rational total_pattern_fraction = id.step_time().fraction();
            if (not id.time_runs_forward()) {
              total_pattern_fraction = 1 - total_pattern_fraction;
            }
            total_pattern_fraction /= abs(pattern_size).fraction();
            return total_pattern_fraction -
                   total_pattern_fraction.numerator() /
                       total_pattern_fraction.denominator();
          };

      const auto time_step =
          local_data->pattern(pattern_fraction_after(local_data->id)) *
          pattern_size;

      local_data->id = local_data->next_id;
      local_data->substep_value = local_data->step_value;
      stepper.add_boundary_delta(make_not_null(&local_data->substep_value),
                                 make_not_null(&local_data->history), time_step,
                                 *local_data->rhs);
      // FIXME dense here?
      if (local_data->id.substep() == 0) {
        local_data->step_value = local_data->substep_value;
      }
      if (local_data->id == final_id) {
        return true;
      }
      const auto next_time_step =
          local_data->pattern(pattern_fraction_after(local_data->next_id)) *
          pattern_size;
      local_data->next_id =
          stepper.next_time_id(local_data->next_id, next_time_step);

      local_data->history.local().insert(local_data->id, stepper.order(),
                                         local_data->substep_value);
      remote_data->messages.emplace_back(local_data->id,
                                         local_data->substep_value);
    }
    return false;
  };

  {
    bool done_x = false;
    bool done_y = false;
    for (;;) {
      done_x = take_steps(&data_x, &data_y);
      if (done_y) {
        REQUIRE(done_x);
        break;
      }
      done_y = take_steps(&data_y, &data_x);
      if (done_x) {
        REQUIRE(done_y);
        break;
      }
    }
  }
  REQUIRE(not dense_output_fraction.has_value());

  CAPTURE(data_x.step_value);
  CAPTURE(data_y.step_value);
  CHECK(data_x.step_value + data_y.step_value == approx(conserved_sum));
  // x + y is conserved, so no need to include y.
  return std::abs(data_x.step_value - expected_x(t_final));
}
}  // namespace

namespace patterns {
struct Lts2to1 {
  static Rational x(const Rational& /*time*/) { return {1, 2}; }
  static Rational y(const Rational& /*time*/) { return {1}; }
};

struct Lts3and1to2 {
  static Rational x(const Rational& time) {
    if (time == 0) {
      return {3, 4};
    } else {
      REQUIRE(time == Rational(3, 4));
      return {1, 4};
    }
  }
  static Rational y(const Rational& /*time*/) { return {1, 2}; }
};
}  // namespace patterns

template <typename Pattern>
void test_convergence(const LtsTimeStepper& stepper,
                      const std::pair<int32_t, int32_t>& step_range,
                      const int32_t stride, const bool output) {
  const auto error = [&](const int32_t intervals) {
    const double forward_error =
        lts_error<Pattern>(stepper, true, intervals, std::nullopt);
    const double backward_error =
        lts_error<Pattern>(stepper, false, intervals, std::nullopt);
    CHECK(forward_error == approx(backward_error));
    return forward_error;
  };
  REQUIRE(TimeStepperTestUtils::convergence_rate(step_range, stride, error,
                                                 output) ==
          approx(stepper.order()).margin(0.4));
}

void test_equal_rate(const LtsTimeStepper& stepper, const size_t order,
                     const size_t number_of_past_steps, const double epsilon,
                     const bool forward) {
  // This does an integral putting the entire derivative into the
  // boundary term.
  const double unused_local_deriv = 4444.;

  auto analytic = [](double t) { return sin(t); };
  auto driver = [](double t) { return cos(t); };
  auto coupling = [=](const double local, const double remote) {
    CHECK(local == unused_local_deriv);
    return remote;
  };

  Approx approx = Approx::custom().epsilon(epsilon);

  const uint64_t num_steps = 100;
  const Slab slab(0.875, 1.);
  const TimeDelta step_size = (forward ? 1 : -1) * slab.duration() / num_steps;

  TimeStepId time_id(forward, 0, forward ? slab.start() : slab.end());
  double y = analytic(time_id.substep_time());
  TimeSteppers::History<double> volume_history{order};
  TimeSteppers::BoundaryHistory<double, double, double> boundary_history{};

  {
    Time history_time = time_id.step_time();
    TimeDelta history_step_size = step_size;
    for (size_t j = 0; j < number_of_past_steps; ++j) {
      ASSERT(history_time.slab() == history_step_size.slab(), "Slab mismatch");
      if ((history_step_size.is_positive() and
           history_time.is_at_slab_start()) or
          (not history_step_size.is_positive() and
           history_time.is_at_slab_end())) {
        const Slab new_slab =
            history_time.slab().advance_towards(-history_step_size);
        history_time = history_time.with_slab(new_slab);
        history_step_size = history_step_size.with_slab(new_slab);
      }
      history_time -= history_step_size;
      const TimeStepId history_id(forward, 0, history_time);
      volume_history.insert_initial(history_id, analytic(history_time.value()),
                                    0.);
      boundary_history.local().insert_initial(history_id, order,
                                              unused_local_deriv);
      boundary_history.remote().insert_initial(history_id, order,
                                               driver(history_time.value()));
    }
  }

  for (uint64_t i = 0; i < num_steps; ++i) {
    for (uint64_t substep = 0;
         substep < stepper.number_of_substeps();
         ++substep) {
      volume_history.insert(time_id, y, 0.);
      boundary_history.local().insert(time_id, order, unused_local_deriv);
      boundary_history.remote().insert(time_id, order,
                                       driver(time_id.substep_time()));

      stepper.update_u(make_not_null(&y), make_not_null(&volume_history),
                       step_size);
      stepper.add_boundary_delta(&y, make_not_null(&boundary_history),
                                 step_size, coupling);
      time_id = stepper.next_time_id(time_id, step_size);
    }
    CHECK(y == approx(analytic(time_id.substep_time())));
  }
  // Make sure history is being cleaned up.  The limit of 20 is
  // arbitrary, but much larger than the order of any integrators we
  // care about and much smaller than the number of time steps in the
  // test.
  CHECK(boundary_history.local().size() < 20);
  CHECK(boundary_history.remote().size() < 20);
}

void test_dense_output(const LtsTimeStepper& stepper) {
  // We only support variable time-step, multistep LTS integration.
  // Any multistep, variable time-step integrator must give the same
  // results from dense output as from just taking a short step
  // because we require dense output to be continuous.  A sufficient
  // test is therefore to run with an LTS pattern and check that the
  // dense output predicts the actual step result.
  const Slab slab(0., 1.);

  // We don't use any meaningful values.  We only care that the dense
  // output gives the same result as normal output.
  // NOLINTNEXTLINE(spectre-mutable)
  auto get_value = [value = 1.]() mutable { return value *= 1.1; };

  const auto coupling = [](const double a, const double b) { return a * b; };

  const auto make_time_id = [](const Time& t) {
    return TimeStepId(true, 0, t);
  };

  TimeSteppers::BoundaryHistory<double, double, double> history{};
  {
    const Slab init_slab = slab.retreat();
    for (size_t i = 0; i < stepper.number_of_past_steps(); ++i) {
      const Time init_time =
          init_slab.end() -
          init_slab.duration() * (i + 1) / stepper.number_of_past_steps();
      history.local().insert_initial(make_time_id(init_time), stepper.order(),
                                     get_value());
      history.remote().insert_initial(make_time_id(init_time), stepper.order(),
                                      get_value());
    }
  }

  std::array<std::deque<TimeDelta>, 2> dt{
      {{slab.duration() / 2, slab.duration() / 4, slab.duration() / 4},
       {slab.duration() / 6, slab.duration() / 6, slab.duration() * 2 / 9,
        slab.duration() * 4 / 9}}};

  Time t = slab.start();
  Time next_check = t + dt[0][0];
  std::array<Time, 2> next{{t, t}};
  for (;;) {
    const size_t side = next[0] <= next[1] ? 0 : 1;

    if (side == 0) {
      history.local().insert(make_time_id(next[0]), stepper.order(),
                             get_value());
    } else {
      history.remote().insert(make_time_id(next[1]), stepper.order(),
                              get_value());
    }

    const TimeDelta this_dt = gsl::at(dt, side).front();
    gsl::at(dt, side).pop_front();

    gsl::at(next, side) += this_dt;

    if (std::min(next[0], next[1]) == next_check) {
      double dense_result = 0.0;
      stepper.boundary_dense_output(&dense_result, history, next_check.value(),
                                    coupling);
      double delta = 0.0;
      stepper.add_boundary_delta(&delta, make_not_null(&history),
                                 next_check - t, coupling);
      CHECK(dense_result == approx(delta));
      if (next_check.is_at_slab_boundary()) {
        break;
      }
      t = next_check;
      next_check += dt[0].front();
    }
  }
}

#define PATTERN(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                               \
  template void test_convergence<PATTERN(data)>(                           \
      const LtsTimeStepper& stepper,                                       \
      const std::pair<int32_t, int32_t>& step_range, const int32_t stride, \
      const bool output);

GENERATE_INSTANTIATIONS(INSTANTIATE, (patterns::Lts2to1, patterns::Lts3and1to2))

#undef INSTANTIATE
#undef PATTERN
}  // namespace TimeStepperTestUtils::lts
