// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <utility>
#include <vector>

#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/Time/TimeSteppers/LtsHelpers.hpp"
#include "Helpers/Time/TimeSteppers/TimeStepperTestUtils.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsMoultonPc.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Literals.hpp"

namespace {
void test_integration() {
  for (size_t order = 2; order < 9; ++order) {
    CAPTURE(order);
    const TimeSteppers::AdamsMoultonPc stepper(order);
    CHECK(stepper.order() == order);
    CHECK(stepper.error_estimate_order() == order - 1);
    CHECK(stepper.number_of_past_steps() == order - 2);
    CHECK(stepper.number_of_substeps() == 2);
    CHECK(stepper.number_of_substeps_for_error() == 2);

    for (size_t start_points = 0;
         start_points <= stepper.number_of_past_steps();
         ++start_points) {
      CAPTURE(start_points);
      const double epsilon = std::max(std::pow(1e-3, start_points + 2), 1e-14);
      TimeStepperTestUtils::integrate_test(stepper, start_points + 2,
                                           start_points, 1., epsilon);
      TimeStepperTestUtils::integrate_test_explicit_time_dependence(
          stepper, start_points + 2, start_points, 1., epsilon);

      const double large_step_epsilon =
          std::clamp(1.0e2 * std::pow(2.0e-2, start_points + 2), 1e-14, 1.0);
      TimeStepperTestUtils::integrate_error_test(
          stepper, start_points + 2, start_points, 1.0, large_step_epsilon, 20,
          1.0e-4);
      TimeStepperTestUtils::integrate_error_test(
          stepper, start_points + 2, start_points, -1.0, large_step_epsilon, 20,
          1.0e-4);

      TimeStepperTestUtils::integrate_variable_test(stepper, start_points + 2,
                                                    start_points, epsilon);
    }
    TimeStepperTestUtils::check_convergence_order(stepper, {10, 30});
    for (size_t history_order = 2; history_order <= order; ++history_order) {
      CAPTURE(history_order);
      TimeStepperTestUtils::check_dense_output(stepper, history_order,
                                               {10, 30});
    }
  }
}

void test_can_change_step_size() {
  const Slab slab(0., 1.);
  const Time start = slab.start();
  const Time mid = slab.start() + slab.duration() / 2;
  const Time end = slab.end();
  const auto can_change = [](const bool time_runs_forward, const Time& first,
                             const Time& second, const Time& now) {
    const TimeSteppers::AdamsMoultonPc stepper(2);
    TimeSteppers::History<double> history(2);
    history.insert(TimeStepId(time_runs_forward, 0, first), 0., 0.);
    history.insert(TimeStepId(time_runs_forward, 2, second), 0., 0.);
    return stepper.can_change_step_size(TimeStepId(time_runs_forward, 4, now),
                                        history);
  };
  CHECK(can_change(true, start, mid, end));
  CHECK_FALSE(can_change(true, start, end, mid));
  CHECK(can_change(true, mid, start, end));
  CHECK_FALSE(can_change(true, mid, end, start));
  CHECK_FALSE(can_change(true, end, start, mid));
  CHECK_FALSE(can_change(true, end, mid, start));
  CHECK(can_change(true, start, mid, mid));
  CHECK_FALSE(can_change(true, start, mid, start));

  CHECK(can_change(false, end, mid, start));
  CHECK_FALSE(can_change(false, end, start, mid));
  CHECK(can_change(false, mid, end, start));
  CHECK_FALSE(can_change(false, mid, start, end));
  CHECK_FALSE(can_change(false, start, end, mid));
  CHECK_FALSE(can_change(false, start, mid, end));
  CHECK(can_change(false, end, mid, mid));
  CHECK_FALSE(can_change(false, end, mid, end));

  {
    const auto time_step = slab.duration() / 2;
    const TimeStepId step_id(true, 1, slab.start() + time_step);
    const TimeSteppers::AdamsMoultonPc stepper(2);
    TimeSteppers::History<double> history(2);
    history.insert(TimeStepId(true, 1, slab.start()), 0., 0.);
    history.insert(step_id, 0., 0.);
    CHECK_FALSE(stepper.can_change_step_size(
        step_id.next_substep(time_step, 1.0), history));
  }
}

void test_equality_and_serialization() {
  TimeSteppers::AdamsMoultonPc am4(4);
  TimeSteppers::AdamsMoultonPc am2(2);
  CHECK(am4 == am4);
  CHECK_FALSE(am4 != am4);
  CHECK(am4 != am2);
  CHECK_FALSE(am4 == am2);

  test_serialization(am4);
  test_serialization_via_base<TimeStepper, TimeSteppers::AdamsMoultonPc>(4_st);
}

void test_creation() {
  const auto created =
      TestHelpers::test_factory_creation<TimeStepper,
                                         TimeSteppers::AdamsMoultonPc>(
          "AdamsMoultonPc:\n"
          "  Order: 3");
  CHECK(created->order() == 3);
}

void test_stability() {
  const auto check_order = [](const size_t order, const double phase) {
    CAPTURE(order);
    TimeStepperTestUtils::stability_test(TimeSteppers::AdamsMoultonPc(order),
                                         phase);
  };

  check_order(2, M_PI);
  check_order(3, 2.504);
  check_order(4, 2.347);
  check_order(5, 2.339);
  check_order(6, 2.368);
  check_order(7, 2.369);
  check_order(8, 2.364);
}

void test_neighbor_data_required() {
  // Test is order-independent
  const TimeSteppers::AdamsMoultonPc stepper(4);
  const Slab slab(0.0, 1.0);
  for (const bool time_runs_forward : {true, false}) {
    const auto start = time_runs_forward ? slab.start() : slab.end();
    const auto duration = (time_runs_forward ? 1 : -1) * slab.duration();
    const auto step_id = [&](const int64_t tenths, const uint64_t substep,
                             const int64_t tenths_step) {
      const auto step_time = start + duration / 10 * tenths;
      if (substep == 0) {
        return TimeStepId(time_runs_forward, 0, step_time);
      } else {
        const auto step_size = duration / 10 * tenths_step;
        return TimeStepId(time_runs_forward, 0, step_time, substep, step_size,
                          (step_time + step_size).value());
      }
    };

    const auto test_ordering = [&](const std::vector<TimeStepId>& ids) {
      for (size_t goal = 0; goal < ids.size(); ++goal) {
        for (size_t data = 0; data < ids.size(); ++data) {
          CHECK(stepper.neighbor_data_required(ids[goal], ids[data]) ==
                (goal > data));
        }
      }
    };

    // GTS
    test_ordering({step_id(0, 0, 1), step_id(0, 1, 1), step_id(1, 0, 1),
                   step_id(1, 1, 1)});

    // LTS
    // 0 1   3
    // 0   2 3
    // step_id(0, 1, 1) and step_id(0, 1, 2) are unsequenced, so we
    // check the sequence with each and then check against each other
    test_ordering({step_id(0, 0, 1), step_id(0, 1, 1), step_id(1, 0, 2),
                   step_id(1, 1, 2), step_id(2, 0, 1), step_id(2, 1, 1),
                   step_id(3, 0, 1)});
    test_ordering({step_id(0, 0, 1), step_id(0, 1, 2), step_id(1, 0, 2),
                   step_id(1, 1, 2), step_id(2, 0, 1), step_id(2, 1, 1),
                   step_id(3, 0, 1)});
    CHECK(
        not stepper.neighbor_data_required(step_id(0, 1, 1), step_id(0, 1, 2)));
    CHECK(
        not stepper.neighbor_data_required(step_id(0, 1, 2), step_id(0, 1, 1)));
  }
}

void test_boundary_gts() {
  for (size_t order = 2; order < 9; ++order) {
    INFO(order);
    const TimeSteppers::AdamsMoultonPc stepper(order);
    for (size_t start_points = 0; start_points < order - 1; ++start_points) {
      INFO(start_points);
      const double epsilon = std::max(std::pow(1e-3, start_points + 2), 1e-14);
      TimeStepperTestUtils::lts::test_equal_rate(stepper, start_points + 2,
                                                 start_points, epsilon, true);
      TimeStepperTestUtils::lts::test_equal_rate(stepper, start_points + 2,
                                                 start_points, epsilon, false);
    }
  }
}

void test_boundary_lts() {
  for (size_t order = 2; order < 9; ++order) {
    INFO(order);
    const TimeSteppers::AdamsMoultonPc stepper(order);

    for (const bool time_runs_forward : {true, false}) {
      const auto error = [&order, &stepper,
                          &time_runs_forward](const int32_t steps) {
        // dx/dt = x y
        // dy/dt = - x y
        //
        // Solution:
        // x = c / [1 + exp(-c (t - d))]
        // y = c - x
        const auto rhs_x = [](const double x, const double y) { return x * y; };
        const auto rhs_y = [](const double x, const double y) {
          return -x * y;
        };

        // Arbitrary except chosen so the error doesn't reach roundoff.
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

        using Hist = TimeSteppers::BoundaryHistory<double, double, double>;
        struct SideData {
          TimeStepId id;
          double value;
          Hist history;
          std::deque<std::pair<TimeStepId, double>> messages;
        };

        const Slab slab =
            time_runs_forward ? Slab(t_init, t_final) : Slab(t_final, t_init);

        SideData data_x;
        SideData data_y;
        {
          // 4.0 to better match step pattern for the evolution.
          const double init_step = std::abs(t_final - t_init) / (4.0 * steps);
          Slab init_slab = slab;
          const auto first_slab = -static_cast<int64_t>(order - 2);
          for (int64_t slab_num = 0; slab_num >= first_slab; --slab_num) {
            const TimeStepId id(
                time_runs_forward, slab_num,
                time_runs_forward ? init_slab.start() : init_slab.end());
            const double x = expected_x(id.step_time().value());
            const double y = conserved_sum - x;
            data_x.history.local().insert_initial(id, order, x);
            data_y.messages.emplace_front(id, x);
            data_y.history.local().insert_initial(id, order, y);
            data_x.messages.emplace_front(id, y);
            init_slab = time_runs_forward
                            ? Slab::with_duration_to_end(
                                  init_slab.start().value(), init_step)
                            : Slab::with_duration_from_start(
                                  init_slab.end().value(), init_step);
          }
        }

        const auto pattern_size =
            (time_runs_forward ? 1 : -1) * slab.duration() / steps;
        const TimeStepId id_final(
            time_runs_forward, 0,
            time_runs_forward ? slab.end() : slab.start());

        data_x.id = TimeStepId(time_runs_forward, 0,
                               time_runs_forward ? slab.start() : slab.end());
        data_y.id = data_x.id;
        data_x.value = expected_x(t_init);
        data_y.value = conserved_sum - data_x.value;

        const auto receive_data =
            [&stepper](const gsl::not_null<SideData*> data, const auto when) {
              while (not data->messages.empty() and
                     stepper.neighbor_data_required(
                         when, data->messages.front().first)) {
                data->history.remote().insert(data->messages.front().first,
                                              stepper.order(),
                                              data->messages.front().second);
                data->messages.pop_front();
              }
            };

        const auto predict = [&receive_data, &stepper](
                                 const gsl::not_null<SideData*> data,
                                 const gsl::not_null<SideData*> other_data,
                                 const TimeDelta& step_size, const auto& rhs) {
          CAPTURE(data->history);
          CAPTURE(step_size);
          data->id = stepper.next_time_id(data->id, step_size);
          receive_data(data, data->id);
          double local_value = data->value;
          stepper.add_boundary_delta(
              &local_value, make_not_null(&data->history), step_size, rhs);
          data->history.local().insert(data->id, stepper.order(), local_value);
          other_data->messages.emplace_back(data->id, local_value);
          // This function is used internally for the corrector.
          if (data->id.substep() == 0) {
            data->value = local_value;
          }
        };

        const auto correct = [&conserved_sum, &predict, &receive_data,
                              &stepper](
                                 const gsl::not_null<SideData*> data,
                                 const gsl::not_null<SideData*> other_data,
                                 const TimeDelta& step_size, const auto& rhs,
                                 const auto& other_rhs) {
          const double time = (data->id.step_time() + step_size).value();
          receive_data(data, time);
          receive_data(other_data, time);
          double dense = data->value;
          double other_dense = other_data->value;
          stepper.boundary_dense_output(&dense, data->history, time, rhs);
          stepper.boundary_dense_output(&other_dense, other_data->history, time,
                                        other_rhs);
          CHECK(dense + other_dense == approx(conserved_sum));
          predict(data, other_data, step_size, rhs);
          CHECK(dense == approx(data->value));
        };

        // Step pattern, repeated `steps` times:
        // 0 1 2 3 4
        // x x     x
        // y     y y
        while (data_x.id != id_final) {
          predict(&data_x, &data_y, pattern_size / 4, rhs_x);
          predict(&data_y, &data_x, pattern_size * 3 / 4, rhs_y);
          correct(&data_x, &data_y, pattern_size / 4, rhs_x, rhs_y);
          predict(&data_x, &data_y, pattern_size * 3 / 4, rhs_x);
          correct(&data_y, &data_x, pattern_size * 3 / 4, rhs_y, rhs_x);
          predict(&data_y, &data_x, pattern_size / 4, rhs_y);
          correct(&data_x, &data_y, pattern_size * 3 / 4, rhs_x, rhs_y);
          correct(&data_y, &data_x, pattern_size / 4, rhs_y, rhs_x);
        }
        // x + y is conserved, so no need to include y.
        return std::abs(data_x.value - expected_x(t_final));
      };
      REQUIRE(TimeStepperTestUtils::convergence_rate(
                  {50 - 4 * order, 130 - 10 * order}, 16, error) ==
              approx(order).margin(0.3));
    }
  }
}

// [[Timeout, 20]]
SPECTRE_TEST_CASE("Unit.Time.TimeSteppers.AdamsMoultonPc", "[Unit][Time]") {
  test_integration();
  test_can_change_step_size();
  test_equality_and_serialization();
  test_creation();
  test_stability();
  test_neighbor_data_required();
  test_boundary_gts();
  test_boundary_lts();
}
}  // namespace
