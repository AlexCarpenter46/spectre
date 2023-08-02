// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Utilities/StdHelpers.hpp"

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <ostream>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "Time/ApproximateTime.hpp"
#include "Time/BoundaryHistory.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsLts.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace {
template <typename T>
void test_apply_coefficients(const T& used_for_size) {
  const Slab slab(1.2, 3.4);

  TimeSteppers::BoundaryHistory<double, double, T> history{};
  history.local().insert(TimeStepId(true, 0, slab.start()), 1, 1.0);
  history.local().insert(TimeStepId(true, 0, slab.end()), 1, 10.0);
  history.remote().insert(TimeStepId(true, 0, slab.start()), 1, 100.0);
  history.remote().insert(TimeStepId(true, 0, slab.end()), 1, 10000.0);

  const TimeSteppers::adams_lts::LtsCoefficients coefficients{
      {history.local()[0], history.remote()[0], 1.0},
      {history.local()[0], history.remote()[1], 2.0},
      {history.local()[1], history.remote()[1], 3.0}};

  const auto coupling = [&used_for_size](const double local,
                                         const double remote) {
    return make_with_value<T>(used_for_size, local * remote);
  };

  auto result = make_with_value<T>(used_for_size, 2.0);
  TimeSteppers::adams_lts::apply_coefficients(
      make_not_null(&result), coefficients, history.evaluator(coupling));

  CHECK(result == make_with_value<T>(used_for_size, 320102.0));
}

class StepChecker {
 public:
  enum class Side { A, B };

  StepChecker(const size_t order) : order_(order) {
    histories_[Side::A];
    histories_[Side::B];
    alt_histories_[Side::A];
    alt_histories_[Side::B];
  }

  void add_step(const Side side, const int step_time) {
    ASSERT(abs(step_time) <= max_time_,
           "Test parameter out of range.  Increase slab size.");

    const TimeStepId id(true, slab_number_, make_time(step_time));
    add_step_impl(side, id);
  }

  void add_substep(const Side side, const int step_time, const int step_size) {
    ASSERT(abs(step_time) <= max_time_,
           "Test parameter out of range.  Increase slab size.");
    ASSERT(abs(step_time + step_size) <= max_time_,
           "Test parameter out of range.  Increase slab size.");

    const TimeStepId id(
        true, slab_number_, make_time(step_time), 1,
        Slab(-max_time_, max_time_).duration() * step_size / (2 * max_time_),
        step_time + step_size);
    add_step_impl(side, id);
  }

  void set_slab_number(const int64_t slab_number) {
    slab_number_ = slab_number;
  }

  struct ExpectedStep {
    int step_time;
    uint64_t substep = 0;
  };
  using ExpectedCoefficients =
      std::map<std::pair<ExpectedStep, ExpectedStep>, double>;

  // This also checks that running backwards with half the step size
  // gives the expected difference.
  ExpectedCoefficients step_coefficients(
      const Side side, const TimeSteppers::adams_lts::StepType step_type,
      const int step_end) const {
    ASSERT(abs(step_end) <= max_time_,
           "Test parameter out of range.  Increase slab size.");

    const auto coefficients =
        coefficients_impl(side, step_type, make_time(step_end));

    if (step_type != TimeSteppers::adams_lts::StepType::Predictor) {
      // Check dense output to the end of the step
      CHECK_ITERABLE_APPROX(
          dense_coefficients(side, step_type, static_cast<double>(step_end)),
          coefficients);
    }

    return coefficients;
  }

  // This also checks that running backwards with half the step size
  // gives the expected difference.
  ExpectedCoefficients dense_coefficients(
      const Side side, const TimeSteppers::adams_lts::StepType step_type,
      const double time) const {
    return coefficients_impl(side, step_type, ApproximateTime{time});
  }

 private:
  static Side other_side(const Side side) {
    return side == Side::A ? Side::B : Side::A;
  }

  static Time make_time(const int time) {
    return Time(Slab(-max_time_, max_time_), {time + max_time_, 2 * max_time_});
  }

  static ExpectedStep expected_index(const TimeStepId& id) {
    return {static_cast<int>(id.step_time().value()), id.substep()};
  }

  void add_step_impl(const Side side, const TimeStepId& id) {
    histories_.at(side).local().insert(id, order_, 0.0);
    histories_.at(other_side(side)).remote().insert(id, order_, 0.0);
    const TimeStepId alt_id = alternate_id(id);
    alt_histories_.at(side).local().insert(alt_id, order_, 0.0);
    alt_histories_.at(other_side(side)).remote().insert(alt_id, order_, 0.0);
  }

  static Time alternate_time(const Time& time) {
    return Time(time.slab(), Time::rational_t{3, 4} - time.fraction() / 2);
  }

  static ApproximateTime alternate_time(const ApproximateTime& time) {
    return ApproximateTime{time.value() / -2.0};
  }

  static TimeStepId alternate_id(const TimeStepId& id) {
    if (id.substep() == 0) {
      return TimeStepId(not id.time_runs_forward(), id.slab_number(),
                        alternate_time(id.step_time()));
    } else {
      return TimeStepId(not id.time_runs_forward(), id.slab_number(),
                        alternate_time(id.step_time()), id.substep(),
                        id.step_size() / -2, id.substep_time() / -2.0);
    }
  }

  static std::map<std::pair<ExpectedStep, ExpectedStep>, double> expected_map(
      const TimeSteppers::adams_lts::LtsCoefficients& coefficients) {
    std::map<std::pair<ExpectedStep, ExpectedStep>, double> coefficients_map{};
    for (const auto& entry : coefficients) {
      const auto insert_success = coefficients_map.insert(
          {{expected_index(get<0>(entry)), expected_index(get<1>(entry))},
           get<2>(entry)});
      if (not insert_success.second) {
        // Duplicate entry.
        CAPTURE(entry);
        CHECK(false);
        insert_success.first->second += get<2>(entry);
      }
    }
    return coefficients_map;
  }

  template <typename TimeType>
  ExpectedCoefficients coefficients_impl(
      const Side side, const TimeSteppers::adams_lts::StepType step_type,
      const TimeType& time) const {
    const auto& history = histories_.at(side);
    const auto coefficients = TimeSteppers::adams_lts::lts_coefficients(
        history.local(), history.remote(), time, step_type);

    {
      const auto& alt_history = alt_histories_.at(side);
      const auto alt_coefficients = TimeSteppers::adams_lts::lts_coefficients(
          alt_history.local(), alt_history.remote(), alternate_time(time),
          step_type);
      REQUIRE(coefficients.size() == alt_coefficients.size());
      for (size_t i = 0; i < coefficients.size(); ++i) {
        const auto& coef = coefficients[i];
        const auto& alt_coef = alt_coefficients[i];
        CHECK(get<0>(alt_coef) == alternate_id(get<0>(coef)));
        CHECK(get<1>(alt_coef) == alternate_id(get<1>(coef)));
        CHECK(get<2>(alt_coef) == approx(get<2>(coef) / -2.0));
      }
    }

    return expected_map(coefficients);
  }

  constexpr static int max_time_ = 16;
  size_t order_;

  int64_t slab_number_{0};
  std::unordered_map<Side,
                     TimeSteppers::BoundaryHistory<double, double, double>>
      histories_;
  std::unordered_map<Side,
                     TimeSteppers::BoundaryHistory<double, double, double>>
      alt_histories_;
};

bool operator<(const StepChecker::ExpectedStep& a,
               const StepChecker::ExpectedStep& b) {
  return a.step_time < b.step_time or
         (a.step_time == b.step_time and a.substep < b.substep);
}

std::ostream& operator<<(std::ostream& s,
                         const StepChecker::ExpectedStep& step) {
  return s << "[" << step.step_time << " " << step.substep << "]";
}

// We call this on the expected values as a check that they make sense.
void verify_conservation(
    const std::vector<StepChecker::ExpectedCoefficients>& coefficients_a,
    const std::vector<StepChecker::ExpectedCoefficients>& coefficients_b) {
  StepChecker::ExpectedCoefficients merged_a{};
  for (const auto& step : coefficients_a) {
    for (const auto& entry : step) {
      const auto insert_success = merged_a.insert(entry);
      if (not insert_success.second) {
        insert_success.first->second += entry.second;
      }
    }
  }

  StepChecker::ExpectedCoefficients merged_b{};
  for (const auto& step : coefficients_b) {
    for (const auto& entry : step) {
      // Swap the sides to make a and b agree.
      const std::pair flipped_index(entry.first.second, entry.first.first);
      const auto insert_success = merged_b.emplace(flipped_index, entry.second);
      if (not insert_success.second) {
        insert_success.first->second += entry.second;
      }
    }
  }

  CHECK_ITERABLE_APPROX(merged_a, merged_b);
}

void test_lts_coefficients() {
  using StepType = TimeSteppers::adams_lts::StepType;
  {
    INFO("AB GTS order 1");
    // 0  1
    // 0  1
    // clang-format off
    const StepChecker::ExpectedCoefficients expected{
        {{{0}, {0}}, 1.0}};
    // clang-format on

    StepChecker checker(1);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 1),
        expected);
  }

  {
    INFO("AB GTS order 3");
    // 0  1  2  3
    // 0  1  2  3
    // clang-format off
    const StepChecker::ExpectedCoefficients expected{
        {{{0}, {0}}, 5.0 / 12.0},
        {{{1}, {1}}, -4.0 / 3.0},
        {{{2}, {2}}, 23.0 / 12.0}};
    // clang-format on

    StepChecker checker(3);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::B, 1);
    checker.add_step(StepChecker::Side::A, 2);
    checker.add_step(StepChecker::Side::B, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 3),
        expected);
  }

  {
    INFO("AB 2:1 order 3");
    // -8          -4           0           4
    //             -4    -2     0     2     4
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_large{
        {{{0}, {2}}, 115.0 / 16.0},
        {{{0}, {0}}, 7.0 / 6.0},
        {{{0}, {-2}}, -11.0 / 16.0},
        {{{-4}, {2}}, -115.0 / 24.0},
        {{{-4}, {-2}}, -11.0 / 8.0},
        {{{-4}, {-4}}, 5.0 / 6.0},
        {{{-8}, {2}}, 23.0 / 16.0},
        {{{-8}, {-2}}, 11.0 / 48.0}};
    const StepChecker::ExpectedCoefficients expected_small_1{
        {{{0}, {0}}, 23.0 / 6.0},
        {{{-2}, {0}}, -1.0},
        {{{-2}, {-4}}, -2.0},
        {{{-4}, {-4}}, 5.0 / 6.0},
        {{{-2}, {-8}}, 1.0 / 3.0}};
    const StepChecker::ExpectedCoefficients expected_small_2{
        {{{2}, {0}}, 115.0 / 16.0},
        {{{0}, {0}}, -8.0 / 3.0},
        {{{-2}, {0}}, 5.0 / 16.0},
        {{{2}, {-4}}, -115.0 / 24.0},
        {{{-2}, {-4}}, 5.0 / 8.0},
        {{{2}, {-8}}, 23.0 / 16.0},
        {{{-2}, {-8}}, -5.0 / 48.0}};
    verify_conservation({expected_large}, {expected_small_1, expected_small_2});
    // clang-format on

    StepChecker checker(3);
    checker.add_step(StepChecker::Side::A, -8);
    checker.add_step(StepChecker::Side::A, -4);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, -4);
    checker.add_step(StepChecker::Side::B, -2);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 2),
        expected_small_1);
    checker.add_step(StepChecker::Side::B, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 4),
        expected_large);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 4),
        expected_small_2);
  }

  {
    INFO("AB LTS -> GTS order 2");
    // -2     0  1
    //    -1  0  1
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_large{
        {{{0}, {0}}, 3.0 / 2.0},
        {{{0}, {-1}}, -1.0 / 4.0},
        {{{-2}, {-1}}, -1.0 / 4.0}};
    const StepChecker::ExpectedCoefficients expected_small{
        {{{0}, {0}}, 3.0 / 2.0},
        {{{-1}, {0}}, -1.0 / 4.0},
        {{{-1}, {-2}}, -1.0 / 4.0}};
    verify_conservation({expected_large}, {expected_small});
    // clang-format on

    StepChecker checker(2);
    checker.add_step(StepChecker::Side::A, -2);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, -1);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 1),
        expected_large);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 1),
        expected_small);
  }

  {
    INFO("AB 3:1 order 2");
    // -3        0        3
    //       -1  0  1  2  3
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_large{
        {{{-3}, {-1}}, -1.0 / 6.0},
        {{{-3}, {1}}, -1.0 / 3.0},
        {{{-3}, {2}}, -1.0},
        {{{0}, {-1}}, -1.0 / 3.0},
        {{{0}, {0}}, 1.0},
        {{{0}, {1}}, 4.0 / 3.0},
        {{{0}, {2}}, 5.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_small_1{
        {{{-1}, {-3}}, -1.0 / 6.0},
        {{{-1}, {0}}, -1.0 / 3.0},
        {{{0}, {0}}, 3.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_small_2{
        {{{0}, {0}}, -1.0 / 2.0},
        {{{1}, {-3}}, -1.0 / 2.0},
        {{{1}, {0}}, 2.0}};
    const StepChecker::ExpectedCoefficients expected_small_3{
        {{{1}, {-3}}, 1.0 / 6.0},
        {{{1}, {0}}, -2.0 / 3.0},
        {{{2}, {-3}}, -1.0},
        {{{2}, {0}}, 5.0 / 2.0}};
    verify_conservation({expected_large},
                        {expected_small_1, expected_small_2, expected_small_3});

    const StepChecker::ExpectedCoefficients expected_dense_large_12{
        {{{-3}, {-1}}, -1.0 / 24.0},
        {{{0}, {-1}}, -1.0 / 12.0},
        {{{0}, {0}}, 5.0 / 8.0}};
    const StepChecker::ExpectedCoefficients expected_dense_small_12{
        {{{-1}, {-3}}, -1.0 / 24.0},
        {{{-1}, {0}}, -1.0 / 12.0},
        {{{0}, {0}}, 5.0 / 8.0}};
    verify_conservation({expected_dense_large_12}, {expected_dense_small_12});
    const StepChecker::ExpectedCoefficients expected_dense_large_32{
        {{{-3}, {-1}}, -1.0 / 6.0},
        {{{-3}, {1}}, -5.0 / 24.0},
        {{{0}, {-1}}, -1.0 / 3.0},
        {{{0}, {0}}, 11.0 / 8.0},
        {{{0}, {1}}, 5.0 / 6.0}};
    const StepChecker::ExpectedCoefficients expected_dense_small_32{
        {{{0}, {0}}, -1.0 / 8.0},
        {{{1}, {-3}}, -5.0 / 24.0},
        {{{1}, {0}}, 5.0 / 6.0}};
    verify_conservation({expected_dense_large_32},
                        {expected_small_1, expected_dense_small_32});
    const StepChecker::ExpectedCoefficients expected_dense_large_52{
        {{{-3}, {-1}}, -1.0 / 6.0},
        {{{-3}, {1}}, -11.0 / 24.0},
        {{{-3}, {2}}, -5.0 / 12.0},
        {{{0}, {-1}}, -1.0 / 3.0},
        {{{0}, {0}}, 1.0},
        {{{0}, {1}}, 11.0 / 6.0},
        {{{0}, {2}}, 25.0 / 24.0}};
    const StepChecker::ExpectedCoefficients expected_dense_small_52{
        {{{1}, {-3}}, 1.0 / 24.0},
        {{{1}, {0}}, -1.0 / 6.0},
        {{{2}, {-3}}, -5.0 / 12.0},
        {{{2}, {0}}, 25.0 / 24.0}};
    verify_conservation(
        {expected_dense_large_52},
        {expected_small_1, expected_small_2, expected_dense_small_52});
    // clang-format on

    StepChecker checker(2);
    checker.add_step(StepChecker::Side::A, -3);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, -1);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::A, StepType::Explicit,
                                   1.0 / 2.0),
        expected_dense_large_12);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::B, StepType::Explicit,
                                   1.0 / 2.0),
        expected_dense_small_12);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 1),
        expected_small_1);
    checker.add_step(StepChecker::Side::B, 1);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::A, StepType::Explicit,
                                   3.0 / 2.0),
        expected_dense_large_32);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::B, StepType::Explicit,
                                   3.0 / 2.0),
        expected_dense_small_32);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 2),
        expected_small_2);
    checker.add_step(StepChecker::Side::B, 2);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::A, StepType::Explicit,
                                   5.0 / 2.0),
        expected_dense_large_52);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::B, StepType::Explicit,
                                   5.0 / 2.0),
        expected_dense_small_52);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 3),
        expected_small_3);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 3),
        expected_large);
  }

  {
    INFO("AB unaligned order 2");
    // 0  1     3  4     6
    // 0     2  3     5  6
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_a_1{
        {{{1}, {2}}, -1.0 / 4.0},
        {{{3}, {2}}, -1.0 / 4.0},
        {{{3}, {3}}, 3.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_a_2{
        {{{3}, {3}}, -1.0 / 2.0},
        {{{3}, {5}}, -3.0 / 2.0},
        {{{4}, {2}}, -3.0 / 2.0},
        {{{4}, {3}}, 11.0 / 4.0},
        {{{4}, {5}}, 11.0 / 4.0}};
    const StepChecker::ExpectedCoefficients expected_b_1{
        {{{2}, {1}}, -1.0 / 4.0},
        {{{2}, {3}}, -1.0 / 4.0},
        {{{2}, {4}}, -3.0 / 2.0},
        {{{3}, {3}}, 1.0},
        {{{3}, {4}}, 3.0}};
    const StepChecker::ExpectedCoefficients expected_b_2{
        {{{3}, {4}}, -1.0 / 4.0},
        {{{5}, {3}}, -3.0 / 2.0},
        {{{5}, {4}}, 11.0 / 4.0}};
    verify_conservation({expected_a_1, expected_a_2},
                        {expected_b_1, expected_b_2});
    // clang-format on

    StepChecker checker(2);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::A, 3);
    checker.add_step(StepChecker::Side::B, 2);
    checker.add_step(StepChecker::Side::B, 3);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 4),
        expected_a_1);
    checker.add_step(StepChecker::Side::A, 4);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 5),
        expected_b_1);
    checker.add_step(StepChecker::Side::B, 5);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 6),
        expected_a_2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Explicit, 6),
        expected_b_2);
  }

  {
    INFO("Predictor-corrector GTS order 2");
    // 0  1
    // 0  1
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_predictor{
        {{{0}, {0}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_corrector{
        {{{0}, {0}}, 1.0 / 2.0},
        {{{0, 1}, {0, 1}}, 1.0 / 2.0}};
    // clang-format on

    StepChecker checker(2);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 1),
        expected_predictor);
    checker.add_substep(StepChecker::Side::A, 0, 1);
    checker.add_substep(StepChecker::Side::B, 0, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 1),
        expected_corrector);
  }

  {
    INFO("Predictor-corrector GTS order 3");
    // 0  1  2
    // 0  1  2
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_predictor{
        {{{0}, {0}}, -1.0 / 2.0},
        {{{1}, {1}}, 3.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_corrector{
        {{{0}, {0}}, -1.0 / 12.0},
        {{{1}, {1}}, 2.0 / 3.0},
        {{{1, 1}, {1, 1}}, 5.0 / 12.0}};
    // clang-format on

    StepChecker checker(3);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::B, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 2),
        expected_predictor);
    checker.add_substep(StepChecker::Side::A, 1, 1);
    checker.add_substep(StepChecker::Side::B, 1, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 2),
        expected_corrector);
  }

  {
    INFO("Predictor-corrector 2:1 order 2");
    // 0     2
    // 0  1  2
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_large_predictor{
        {{{0}, {0}}, 2.0}};
    const StepChecker::ExpectedCoefficients expected_large_corrector{
        {{{0}, {0}}, 1.0 / 2.0},
        {{{0}, {0, 1}}, 1.0 / 4.0},
        {{{0, 1}, {0, 1}}, 1.0 / 4.0},
        {{{0}, {1}}, 1.0 / 4.0},
        {{{0, 1}, {1}}, 1.0 / 4.0},
        {{{0, 1}, {1, 1}}, 1.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_small_predictor_1{
        {{{0}, {0}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_small_corrector_1{
        {{{0}, {0}}, 1.0 / 2.0},
        {{{0, 1}, {0}}, 1.0 / 4.0},
        {{{0, 1}, {0, 1}}, 1.0 / 4.0}};
    const StepChecker::ExpectedCoefficients expected_small_predictor_2{
        {{{1}, {0, 1}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_small_corrector_2{
        {{{1}, {0}}, 1.0 / 4.0},
        {{{1}, {0, 1}}, 1.0 / 4.0},
        {{{1, 1}, {0, 1}}, 1.0 / 2.0}};
    verify_conservation(
        {expected_large_corrector},
        {expected_small_corrector_1, expected_small_corrector_2});

    const StepChecker::ExpectedCoefficients expected_dense_large_12{
        {{{0}, {0}}, 3.0 / 8.0},
        {{{0}, {0, 1}}, 1.0 / 16.0},
        {{{0, 1}, {0, 1}}, 1.0 / 16.0}};
    const StepChecker::ExpectedCoefficients expected_dense_small_12{
        {{{0}, {0}}, 3.0 / 8.0},
        {{{0, 1}, {0}}, 1.0 / 16.0},
        {{{0, 1}, {0, 1}}, 1.0 / 16.0}};
    verify_conservation({expected_dense_large_12}, {expected_dense_small_12});
    const StepChecker::ExpectedCoefficients expected_dense_large_32{
        {{{0}, {0}}, 1.0 / 2.0},
        {{{0}, {0, 1}}, 1.0 / 4.0},
        {{{0}, {1}}, 3.0 / 16.0},
        {{{0, 1}, {0, 1}}, 1.0 / 4.0},
        {{{0, 1}, {1}}, 3.0 / 16.0},
        {{{0, 1}, {1, 1}}, 1.0 / 8.0}};
    const StepChecker::ExpectedCoefficients expected_dense_small_32{
        {{{1}, {0}}, 3.0 / 16.0},
        {{{1}, {0, 1}}, 3.0 / 16.0},
        {{{1, 1}, {0, 1}}, 1.0 / 8.0}};
    verify_conservation({expected_dense_large_32},
                        {expected_small_corrector_1, expected_dense_small_32});
    // clang-format on

    StepChecker checker(2);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 2),
        expected_large_predictor);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 1),
        expected_small_predictor_1);
    checker.add_substep(StepChecker::Side::A, 0, 2);
    checker.add_substep(StepChecker::Side::B, 0, 1);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::A, StepType::Corrector,
                                   1.0 / 2.0),
        expected_dense_large_12);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::B, StepType::Corrector,
                                   1.0 / 2.0),
        expected_dense_small_12);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 1),
        expected_small_corrector_1);
    checker.add_step(StepChecker::Side::B, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 2),
        expected_small_predictor_2);
    checker.add_substep(StepChecker::Side::B, 1, 1);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::A, StepType::Corrector,
                                   3.0 / 2.0),
        expected_dense_large_32);
    CHECK_ITERABLE_APPROX(
        checker.dense_coefficients(StepChecker::Side::B, StepType::Corrector,
                                   3.0 / 2.0),
        expected_dense_small_32);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 2),
        expected_small_corrector_2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 2),
        expected_large_corrector);
  }

  {
    INFO("Predictor-corrector 2:1 order 3");
    // 0     2     4
    //    1  2  3  4
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_large_predictor{
        {{{0}, {1}}, -1.0},
        {{{2}, {1}}, -1.0},
        {{{2}, {2}}, 4.0}};
    const StepChecker::ExpectedCoefficients expected_large_corrector{
        {{{0}, {1}}, -1.0 / 32.0},
        {{{0}, {2, 1}}, -5.0 / 96.0},
        {{{0}, {3}}, -1.0 / 12.0},
        {{{2}, {1}}, -1.0 / 16.0},
        {{{2}, {2}}, 7.0 / 12.0},
        {{{2}, {2, 1}}, 5.0 / 16.0},
        {{{2}, {3}}, 1.0 / 2.0},
        {{{2, 1}, {1}}, 1.0 / 96.0},
        {{{2, 1}, {2, 1}}, 5.0 / 32.0},
        {{{2, 1}, {3}}, 1.0 / 4.0},
        {{{2, 1}, {3, 1}}, 5.0 / 12.0}};
    const StepChecker::ExpectedCoefficients expected_small_predictor_1{
        {{{1}, {0}}, -1.0 / 4.0},
        {{{1}, {2}}, -1.0 / 4.0},
        {{{2}, {2}}, 3.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_small_corrector_1{
        {{{1}, {0}}, -1.0 / 32.0},
        {{{1}, {2}}, -1.0 / 16.0},
        {{{1}, {2, 1}}, 1.0 / 96.0},
        {{{2}, {2}}, 2.0 / 3.0},
        {{{2, 1}, {0}}, -5.0 / 96.0},
        {{{2, 1}, {2}}, 5.0 / 16.0},
        {{{2, 1}, {2, 1}}, 5.0 / 32.0}};
    const StepChecker::ExpectedCoefficients expected_small_predictor_2{
        {{{2}, {2}}, -1.0 / 2.0},
        {{{3}, {2}}, 3.0 / 4.0},
        {{{3}, {2, 1}}, 3.0 / 4.0}};
    const StepChecker::ExpectedCoefficients expected_small_corrector_2{
        {{{2}, {2}}, -1.0 / 12.0},
        {{{3}, {0}}, -1.0 / 12.0},
        {{{3}, {2}}, 1.0 / 2.0},
        {{{3}, {2, 1}}, 1.0 / 4.0},
        {{{3, 1}, {2, 1}}, 5.0 / 12.0}};
    verify_conservation(
        {expected_large_corrector},
        {expected_small_corrector_1, expected_small_corrector_2});
    // clang-format on

    StepChecker checker(3);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 1);
    checker.add_step(StepChecker::Side::A, 2);
    checker.add_step(StepChecker::Side::B, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 4),
        expected_large_predictor);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 3),
        expected_small_predictor_1);
    checker.add_substep(StepChecker::Side::A, 2, 2);
    checker.add_substep(StepChecker::Side::B, 2, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 3),
        expected_small_corrector_1);
    checker.add_step(StepChecker::Side::B, 3);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 4),
        expected_small_predictor_2);
    checker.add_substep(StepChecker::Side::B, 3, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 4),
        expected_small_corrector_2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 4),
        expected_large_corrector);
  }

  {
    INFO("Predictor-corrector 3:1 order 2");
    // 0        3
    // 0  1  2  3
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_large_predictor{
        {{{0}, {0}}, 3.0}};
    const StepChecker::ExpectedCoefficients expected_large_corrector{
        {{{0}, {0}}, 1.0 / 2.0},
        {{{0}, {0, 1}}, 1.0 / 3.0},
        {{{0}, {1}}, 1.0 / 3.0},
        {{{0}, {1, 1}}, 1.0 / 6.0},
        {{{0}, {2}}, 1.0 / 6.0},
        {{{0, 1}, {0, 1}}, 1.0 / 6.0},
        {{{0, 1}, {1}}, 1.0 / 6.0},
        {{{0, 1}, {1, 1}}, 1.0 / 3.0},
        {{{0, 1}, {2}}, 1.0 / 3.0},
        {{{0, 1}, {2, 1}}, 1.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_small_predictor_1{
        {{{0}, {0}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_small_corrector_1{
        {{{0}, {0}}, 1.0 / 2.0},
        {{{0, 1}, {0}}, 1.0 / 3.0},
        {{{0, 1}, {0, 1}}, 1.0 / 6.0}};
    const StepChecker::ExpectedCoefficients expected_small_predictor_2{
        {{{1}, {0, 1}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_small_corrector_2{
        {{{1}, {0}}, 1.0 / 3.0},
        {{{1}, {0, 1}}, 1.0 / 6.0},
        {{{1, 1}, {0}}, 1.0 / 6.0},
        {{{1, 1}, {0, 1}}, 1.0 / 3.0}};
    const StepChecker::ExpectedCoefficients expected_small_predictor_3{
        {{{2}, {0, 1}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_small_corrector_3{
        {{{2}, {0}}, 1.0 / 6.0},
        {{{2}, {0, 1}}, 1.0 / 3.0},
        {{{2, 1}, {0, 1}}, 1.0 / 2.0}};
    verify_conservation({expected_large_corrector},
                        {expected_small_corrector_1, expected_small_corrector_2,
                         expected_small_corrector_3});
    // clang-format on

    StepChecker checker(2);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 3),
        expected_large_predictor);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 1),
        expected_small_predictor_1);
    checker.add_substep(StepChecker::Side::A, 0, 3);
    checker.add_substep(StepChecker::Side::B, 0, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 1),
        expected_small_corrector_1);
    checker.add_step(StepChecker::Side::B, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 2),
        expected_small_predictor_2);
    checker.add_substep(StepChecker::Side::B, 1, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 2),
        expected_small_corrector_2);
    checker.add_step(StepChecker::Side::B, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 3),
        expected_small_predictor_3);
    checker.add_substep(StepChecker::Side::B, 2, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 3),
        expected_small_corrector_3);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 3),
        expected_large_corrector);
  }

  {
    INFO("Predictor-corrector unaligned order 2");
    //          3  4     6
    //          3     5  6
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_a_predictor_1{
        {{{3}, {3}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_a_corrector_1{
        {{{3}, {3}}, 1.0 / 2.0},
        {{{3, 1}, {3}}, 1.0 / 4.0},
        {{{3, 1}, {3, 1}}, 1.0 / 4.0}};
    const StepChecker::ExpectedCoefficients expected_a_predictor_2{
        {{{4}, {3, 1}}, 2.0}};
    const StepChecker::ExpectedCoefficients expected_a_corrector_2{
        {{{4}, {3}}, 1.0 / 4.0},
        {{{4}, {3, 1}}, 1.0 / 2.0},
        {{{4}, {5}}, 1.0 / 4.0},
        {{{4, 1}, {3, 1}}, 1.0 / 4.0},
        {{{4, 1}, {5}}, 1.0 / 4.0},
        {{{4, 1}, {5, 1}}, 1.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_b_predictor_1{
        {{{3}, {3}}, 2.0}};
    const StepChecker::ExpectedCoefficients expected_b_corrector_1{
        {{{3}, {3}}, 1.0 / 2.0},
        {{{3}, {3, 1}}, 1.0 / 4.0},
        {{{3}, {4}}, 1.0 / 4.0},
        {{{3, 1}, {3, 1}}, 1.0 / 4.0},
        {{{3, 1}, {4}}, 1.0 / 2.0},
        {{{3, 1}, {4, 1}}, 1.0 / 4.0}};
    const StepChecker::ExpectedCoefficients expected_b_predictor_2{
        {{{5}, {4, 1}}, 1.0}};
    const StepChecker::ExpectedCoefficients expected_b_corrector_2{
        {{{5}, {4}}, 1.0 / 4.0},
        {{{5}, {4, 1}}, 1.0 / 4.0},
        {{{5, 1}, {4, 1}}, 1.0 / 2.0}};
    verify_conservation({expected_a_corrector_1, expected_a_corrector_2},
                        {expected_b_corrector_1, expected_b_corrector_2});
    // clang-format on

    StepChecker checker(2);
    checker.add_step(StepChecker::Side::A, 3);
    checker.add_step(StepChecker::Side::B, 3);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 4),
        expected_a_predictor_1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 5),
        expected_b_predictor_1);
    checker.add_substep(StepChecker::Side::A, 3, 1);
    checker.add_substep(StepChecker::Side::B, 3, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 4),
        expected_a_corrector_1);
    checker.add_step(StepChecker::Side::A, 4);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 6),
        expected_a_predictor_2);
    checker.add_substep(StepChecker::Side::A, 4, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 5),
        expected_b_corrector_1);
    checker.add_step(StepChecker::Side::B, 5);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 6),
        expected_b_predictor_2);
    checker.add_substep(StepChecker::Side::B, 5, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 6),
        expected_a_corrector_2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 6),
        expected_b_corrector_2);
  }

  {
    INFO("Predictor-corrector unaligned order 3");
    //    1     3  4     6
    //       2  3     5  6
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_a_predictor_1{
        {{{1}, {2}}, -1.0 / 4.0},
        {{{3}, {2}}, -1.0 / 4.0},
        {{{3}, {3}}, 3.0 / 2.0}};
    const StepChecker::ExpectedCoefficients expected_a_corrector_1{
        {{{1}, {2}}, -1.0 / 36.0},
        {{{3}, {2}}, -1.0 / 12.0},
        {{{3}, {3}}, 2.0 / 3.0},
        {{{3, 1}, {2}}, -1.0 / 9.0},
        {{{3, 1}, {3}}, 5.0 / 12.0},
        {{{3, 1}, {3, 1}}, 5.0 / 36.0}};
    const StepChecker::ExpectedCoefficients expected_a_predictor_2{
        {{{3}, {3}}, -2.0},
        {{{4}, {3}}, 2.0},
        {{{4}, {3, 1}}, 2.0}};
    const StepChecker::ExpectedCoefficients expected_a_corrector_2{
        {{{3}, {3}}, -1.0 / 12.0},
        {{{3}, {3, 1}}, -5.0 / 36.0},
        {{{3}, {5}}, -2.0 / 9.0},
        {{{4}, {2}}, -2.0 / 9.0},
        {{{4}, {3}}, 23.0 / 36.0},
        {{{4}, {3, 1}}, 23.0 / 36.0},
        {{{4}, {5}}, 7.0 / 12.0},
        {{{4}, {5, 1}}, 1.0 / 36.0},
        {{{4, 1}, {3, 1}}, 5.0 / 36.0},
        {{{4, 1}, {5}}, 2.0 / 9.0},
        {{{4, 1}, {5, 1}}, 5.0 / 12.0}};
    const StepChecker::ExpectedCoefficients expected_b_predictor_1{
        {{{2}, {1}}, -1.0},
        {{{2}, {3}}, -1.0},
        {{{3}, {3}}, 4.0}};
    const StepChecker::ExpectedCoefficients expected_b_corrector_1{
        {{{2}, {1}}, -1.0 / 36.0},
        {{{2}, {3}}, -1.0 / 12.0},
        {{{2}, {3, 1}}, -1.0 / 9.0},
        {{{2}, {4}}, -2.0 / 9.0},
        {{{3}, {3}}, 7.0 / 12.0},
        {{{3}, {3, 1}}, 5.0 / 12.0},
        {{{3}, {4}}, 2.0 / 3.0},
        {{{3, 1}, {3}}, -5.0 / 36.0},
        {{{3, 1}, {3, 1}}, 5.0 / 36.0},
        {{{3, 1}, {4}}, 23.0 / 36.0},
        {{{3, 1}, {4, 1}}, 5.0 / 36.0}};
    const StepChecker::ExpectedCoefficients expected_b_predictor_2{
        {{{3}, {4}}, -1.0 / 4.0},
        {{{5}, {4}}, 1.0 / 2.0},
        {{{5}, {4, 1}}, 3.0 / 4.0}};
    const StepChecker::ExpectedCoefficients expected_b_corrector_2{
        {{{3}, {4}}, -1.0 / 36.0},
        {{{5}, {3}}, -2.0 / 9.0},
        {{{5}, {4}}, 7.0 / 12.0},
        {{{5}, {4, 1}}, 2.0 / 9.0},
        {{{5, 1}, {4}}, 1.0 / 36.0},
        {{{5, 1}, {4, 1}}, 5.0 / 12.0}};
    verify_conservation({expected_a_corrector_1, expected_a_corrector_2},
                        {expected_b_corrector_1, expected_b_corrector_2});
    // clang-format on

    StepChecker checker(3);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::B, 2);
    checker.add_step(StepChecker::Side::A, 3);
    checker.add_step(StepChecker::Side::B, 3);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 4),
        expected_a_predictor_1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 5),
        expected_b_predictor_1);
    checker.add_substep(StepChecker::Side::A, 3, 1);
    checker.add_substep(StepChecker::Side::B, 3, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 4),
        expected_a_corrector_1);
    checker.add_step(StepChecker::Side::A, 4);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 6),
        expected_a_predictor_2);
    checker.add_substep(StepChecker::Side::A, 4, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 5),
        expected_b_corrector_1);
    checker.add_step(StepChecker::Side::B, 5);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Predictor, 6),
        expected_b_predictor_2);
    checker.add_substep(StepChecker::Side::B, 5, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 6),
        expected_a_corrector_2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::B, StepType::Corrector, 6),
        expected_b_corrector_2);
  }

  {
    INFO("AB self-start order 3");
    // 0 1 2 -- 0 1 2
    // 0 1 2 -- 0 1 2
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_1{
        {{{0}, {0}}, 5.0 / 12.0},
        {{{1}, {1}}, 2.0 / 3.0},
        {{{2}, {2}}, -1.0 / 12.0}};
    const StepChecker::ExpectedCoefficients expected_2{
        {{{0}, {0}}, -1.0 / 12.0},
        {{{1}, {1}}, 2.0 / 3.0},
        {{{2}, {2}}, 5.0 / 12.0}};
    // clang-format on

    StepChecker checker(3);
    checker.set_slab_number(-1);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::B, 1);
    checker.add_step(StepChecker::Side::A, 2);
    checker.add_step(StepChecker::Side::B, 2);
    checker.set_slab_number(0);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 1),
        expected_1);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::B, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Explicit, 2),
        expected_2);
  }

  {
    INFO("Predictor-corrector self-start order 4");
    // 0 1 2 3 -- 0 1 2 3
    // 0 1 2 3 -- 0 1 2 3
    // clang-format off
    const StepChecker::ExpectedCoefficients expected_predictor_1{
        {{{0}, {0}}, 23.0 / 36.0},
        {{{2}, {2}}, 7.0 / 12.0},
        {{{3}, {3}}, -2.0 / 9.0}};
    const StepChecker::ExpectedCoefficients expected_corrector_1{
        {{{0}, {0}}, 3.0 / 8.0},
        {{{0, 1}, {0, 1}}, 19.0 / 24.0},
        {{{2}, {2}}, -5.0 / 24.0},
        {{{3}, {3}}, 1.0 / 24.0}};
    const StepChecker::ExpectedCoefficients expected_predictor_2{
        {{{0}, {0}}, -2.0 / 9.0},
        {{{1}, {1}}, 13.0 / 12.0},
        {{{3}, {3}}, 5.0 / 36.0}};
    const StepChecker::ExpectedCoefficients expected_corrector_2{
        {{{0}, {0}}, -1.0 / 24.0},
        {{{1}, {1}}, 13.0 / 24.0},
        {{{1, 1}, {1, 1}}, 13.0 / 24.0},
        {{{3}, {3}}, -1.0 / 24.0}};
    const StepChecker::ExpectedCoefficients expected_predictor_3{
        {{{0}, {0}}, 5.0 / 12.0},
        {{{1}, {1}}, -4.0 / 3.0},
        {{{2}, {2}}, 23.0 / 12.0}};
    const StepChecker::ExpectedCoefficients expected_corrector_3{
        {{{0}, {0}}, 1.0 / 24.0},
        {{{1}, {1}}, -5.0 / 24.0},
        {{{2}, {2}}, 19.0 / 24.0},
        {{{2, 1}, {2, 1}}, 3.0 / 8.0}};
    // clang-format on

    StepChecker checker(4);
    checker.set_slab_number(-1);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::B, 1);
    checker.add_step(StepChecker::Side::A, 2);
    checker.add_step(StepChecker::Side::B, 2);
    checker.add_step(StepChecker::Side::A, 3);
    checker.add_step(StepChecker::Side::B, 3);
    checker.set_slab_number(0);
    checker.add_step(StepChecker::Side::A, 0);
    checker.add_step(StepChecker::Side::B, 0);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 1),
        expected_predictor_1);
    checker.add_substep(StepChecker::Side::A, 0, 1);
    checker.add_substep(StepChecker::Side::B, 0, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 1),
        expected_corrector_1);
    checker.add_step(StepChecker::Side::A, 1);
    checker.add_step(StepChecker::Side::B, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 2),
        expected_predictor_2);
    checker.add_substep(StepChecker::Side::A, 1, 1);
    checker.add_substep(StepChecker::Side::B, 1, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 2),
        expected_corrector_2);
    checker.add_step(StepChecker::Side::A, 2);
    checker.add_step(StepChecker::Side::B, 2);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Predictor, 3),
        expected_predictor_3);
    checker.add_substep(StepChecker::Side::A, 2, 1);
    checker.add_substep(StepChecker::Side::B, 2, 1);
    CHECK_ITERABLE_APPROX(
        checker.step_coefficients(StepChecker::Side::A, StepType::Corrector, 3),
        expected_corrector_3);
  }
}

void test_clean_boundary_history() {
  for (const bool time_runs_forward : {true, false}) {
    CAPTURE(time_runs_forward);
    const auto make_id = [&time_runs_forward](const int32_t step,
                                              const uint64_t substep,
                                              const int32_t step_scale) {
      const Slab slab(0.0, 1.0);
      const auto base_step =
          (time_runs_forward ? 1 : -1) * slab.duration() / 20;
      const auto time = time_runs_forward ? slab.start() + step * base_step
                                          : slab.end() + step * base_step;
      const auto step_size = base_step * step_scale;
      return TimeStepId(time_runs_forward, 0, time, substep, step_size,
                        (time + substep * step_size).value());
    };

    // GTS Euler
    {
      TimeSteppers::BoundaryHistory<double, double, double> hist{};
      // The order from the history should be ignored.
      hist.local().insert(make_id(0, 0, 1), 5, 0.0);
      hist.remote().insert(make_id(0, 0, 1), 5, 0.0);
      TimeSteppers::adams_lts::clean_boundary_history(hist.local(),
                                                      hist.remote(), 1);
      CHECK(hist.local().size() == 1);
      CHECK(hist.remote().size() == 1);

      hist.local().insert(make_id(1, 0, 1), 5, 0.0);
      hist.remote().insert(make_id(1, 0, 1), 5, 0.0);
      TimeSteppers::adams_lts::clean_boundary_history(hist.local(),
                                                      hist.remote(), 1);
      CHECK(hist.local().size() == 1);
      CHECK(hist.remote().size() == 1);
    }

    // GTS
    {
      TimeSteppers::BoundaryHistory<double, double, double> hist{};
      for (int32_t i = 0; i <= 8; ++i) {
        hist.local().insert(make_id(i, 0, 1), 5, 0.0);
        hist.local().insert(make_id(i, 1, 1), 5, 0.0);
        hist.remote().insert(make_id(i, 0, 1), 5, 0.0);
        hist.remote().insert(make_id(i, 1, 1), 5, 0.0);
      }
      hist.local().clear_substeps(hist.local().size() - 1);
      hist.remote().clear_substeps(hist.remote().size() - 1);
      TimeSteppers::adams_lts::clean_boundary_history(hist.local(),
                                                      hist.remote(), 3);
      CHECK(hist.local().size() == 3);
      CHECK(hist.remote().size() == 3);
      CHECK(hist.local().number_of_substeps(0) == 1);
      CHECK(hist.local().number_of_substeps(1) == 1);
      CHECK(hist.local().number_of_substeps(2) == 1);
      CHECK(hist.remote().number_of_substeps(0) == 1);
      CHECK(hist.remote().number_of_substeps(1) == 1);
      CHECK(hist.remote().number_of_substeps(2) == 1);
    }

    // LTS local densest
    {
      TimeSteppers::BoundaryHistory<double, double, double> hist{};
      for (int32_t i = 0; i <= 8; ++i) {
        hist.local().insert(make_id(i, 0, 1), 5, 0.0);
        hist.local().insert(make_id(i, 1, 1), 5, 0.0);
        if (i % 3 == 0) {
          hist.remote().insert(make_id(i, 0, 3), 5, 0.0);
          hist.remote().insert(make_id(i, 1, 3), 5, 0.0);
        }
      }
      hist.local().clear_substeps(hist.local().size() - 1);
      TimeSteppers::adams_lts::clean_boundary_history(hist.local(),
                                                      hist.remote(), 2);
      CHECK(hist.local().size() == 2);
      CHECK(hist.remote().size() == 2);
      CHECK(hist.local().number_of_substeps(0) == 1);
      CHECK(hist.local().number_of_substeps(1) == 1);
      CHECK(hist.remote().number_of_substeps(0) == 1);
      CHECK(hist.remote().number_of_substeps(1) == 2);
    }

    // LTS remote densest
    {
      TimeSteppers::BoundaryHistory<double, double, double> hist{};
      for (int32_t i = 0; i <= 8; ++i) {
        if (i % 3 == 0) {
          hist.local().insert(make_id(i, 0, 3), 5, 0.0);
          hist.local().insert(make_id(i, 1, 3), 5, 0.0);
        }
        hist.remote().insert(make_id(i, 0, 1), 5, 0.0);
        hist.remote().insert(make_id(i, 1, 1), 5, 0.0);
      }
      hist.local().clear_substeps(hist.local().size() - 1);
      hist.remote().clear_substeps(hist.remote().size() - 1);
      TimeSteppers::adams_lts::clean_boundary_history(hist.local(),
                                                      hist.remote(), 2);
      CHECK(hist.local().size() == 2);
      // The remote side needs more data for the earlier small steps.
      CHECK(hist.remote().size() == 4);
      CHECK(hist.local().number_of_substeps(0) == 1);
      CHECK(hist.local().number_of_substeps(1) == 1);
      CHECK(hist.remote().number_of_substeps(0) == 1);
      CHECK(hist.remote().number_of_substeps(1) == 2);
      CHECK(hist.remote().number_of_substeps(2) == 2);
      CHECK(hist.remote().number_of_substeps(3) == 1);
    }

    // LTS interleaved - predictor is in interior of remote step
    {
      TimeSteppers::BoundaryHistory<double, double, double> hist{};
      hist.local().insert(make_id(0, 0, 2), 5, 0.0);
      hist.remote().insert(make_id(0, 0, 1), 5, 0.0);
      hist.remote().insert(make_id(0, 1, 1), 5, 0.0);
      hist.remote().insert(make_id(1, 0, 2), 5, 0.0);
      hist.remote().insert(make_id(1, 1, 2), 5, 0.0);
      TimeSteppers::adams_lts::clean_boundary_history(hist.local(),
                                                      hist.remote(), 1);
      CHECK(hist.local().size() == 1);
      CHECK(hist.remote().size() == 2);
      CHECK(hist.local().number_of_substeps(0) == 1);
      CHECK(hist.remote().number_of_substeps(0) == 2);
      CHECK(hist.remote().number_of_substeps(1) == 2);
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Time.TimeSteppers.AdamsLts", "[Unit][Time]") {
  test_apply_coefficients(0.0);
  test_apply_coefficients(DataVector(5, 0.0));
  test_lts_coefficients();
  test_clean_boundary_history();
}
