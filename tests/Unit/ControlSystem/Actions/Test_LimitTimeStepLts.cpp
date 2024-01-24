// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ControlSystem/Actions/LimitTimeStepLts.hpp"
#include "ControlSystem/FutureMeasurements.hpp"
#include "ControlSystem/Tags/FutureMeasurements.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/AdaptiveSteppingDiagnostics.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsBashforth.hpp"
#include "Time/TimeSteppers/AdamsMoultonPc.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

// FIXME can we reduce copy-paste from the GTS version?
namespace {
struct LabelA {};
struct LabelB {};
struct LabelC {};

using MeasurementA = control_system::TestHelpers::Measurement<LabelA>;
using MeasurementB = control_system::TestHelpers::Measurement<LabelB>;

using systemsA =
    tmpl::list<control_system::TestHelpers::System<0, LabelA, MeasurementA>>;
using systemsB =
    tmpl::list<control_system::TestHelpers::System<0, LabelB, MeasurementB>,
               control_system::TestHelpers::System<0, LabelC, MeasurementB>>;

using control_systems = tmpl::append<systemsA, systemsB>;

struct Var : db::SimpleTag {
  using type = double;
};

struct FunctionsOfTimeTag : domain::Tags::FunctionsOfTime, db::SimpleTag {
  using type = std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;
};

template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = ElementId<1>;
  using mutable_global_cache_tags =
      tmpl::list<control_system::Tags::MeasurementTimescales,
                 FunctionsOfTimeTag>;
  // FIXME is everything needed?
  using simple_tags = db::AddSimpleTags<
      domain::Tags::Domain<1>,
      control_system::Tags::FutureMeasurements<systemsA>,
      control_system::Tags::FutureMeasurements<systemsB>,
      Tags::ConcreteTimeStepper<LtsTimeStepper>, Tags::TimeStepId,
      Tags::Next<Tags::TimeStepId>, Tags::TimeStep, Tags::Next<Tags::TimeStep>,
      Tags::AdaptiveSteppingDiagnostics, Tags::HistoryEvolvedVariables<Var>>;
  using compute_tags = time_stepper_ref_tags<LtsTimeStepper>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 simple_tags, compute_tags>>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<
              control_system::Actions::LimitTimeStepLts<control_systems>>>>;
};

struct Metavariables {
  static constexpr size_t volume_dim = 1;
  static constexpr bool local_time_stepping = true;
  using component_list = tmpl::list<Component<Metavariables>>;
};

void test_impl(
    const std::string& test_label, const bool in_measurement_region,
    const Rational initial_time, const Rational initial_step_end,
    const std::optional<Rational>& expected_step_end,
    const std::vector<std::pair<double, double>>& measurement_updatesA,
    const std::vector<std::pair<double, double>>& fot_updatesA,
    const std::vector<std::pair<double, double>>& measurement_updatesBC,
    const std::vector<std::pair<double, double>>& fot_updatesB,
    const std::vector<std::pair<double, double>>& fot_updatesC,
    std::unique_ptr<LtsTimeStepper> stepper,
    TimeSteppers::History<Var::type> history) {
  INFO(test_label);
  ASSERT(measurement_updatesA.size() > 1, "Bad argument");
  ASSERT(measurement_updatesBC.size() > 1, "Bad argument");
  ASSERT(fot_updatesA.size() > 1, "Bad argument");
  ASSERT(fot_updatesB.size() > 1, "Bad argument");
  ASSERT(fot_updatesC.size() > 1, "Bad argument");

  std::vector<Block<1>> blocks{};
  {
    Block<1> block(
        domain::make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
            domain::CoordinateMaps::Identity<1>{}),
        1, {});
    if (in_measurement_region) {
      block.inject_time_dependent_map(
          domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
              domain::CoordinateMaps::Identity<1>{}),
          domain::make_coordinate_map_base<Frame::Grid, Frame::Distorted>(
              domain::CoordinateMaps::Identity<1>{}),
          domain::make_coordinate_map_base<Frame::Distorted, Frame::Inertial>(
              domain::CoordinateMaps::Identity<1>{}));
    }
    blocks.push_back(std::move(block));
  }
  Domain<1> domain(std::move(blocks));
  const ElementId<1> element_id(0);

  const Slab slab(0.0, 1.0);
  const TimeStepId initial_id(true, 0,
                              slab.start() + initial_time * slab.duration());

  const size_t measurements_per_update = 3;

  const auto setup_fot =
      [](const std::vector<std::pair<double, double>>& updates) {
        auto fot =
            std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
                updates.front().first,
                std::array{DataVector{updates.front().second}},
                updates[1].first);
        for (size_t i = 1; i < updates.size() - 1; ++i) {
          fot->update(updates[i].first, DataVector{updates[i].second},
                      updates[i + 1].first);
        }
        return fot;
      };

  control_system::Tags::MeasurementTimescales::type timescales{};
  timescales["LabelA"] = setup_fot(measurement_updatesA);
  timescales["LabelBLabelC"] = setup_fot(measurement_updatesBC);
  FunctionsOfTimeTag::type functions_of_time{};
  functions_of_time["LabelA"] = setup_fot(fot_updatesA);
  functions_of_time["LabelB"] = setup_fot(fot_updatesB);
  functions_of_time["LabelC"] = setup_fot(fot_updatesC);

  const auto setup_measurements =
      [&initial_id](const domain::FunctionsOfTime::FunctionOfTime& timescale) {
        if (timescale.func(0.0)[0][0] ==
            std::numeric_limits<double>::infinity()) {
          return control_system::FutureMeasurements(
              1, std::numeric_limits<double>::infinity());
        }

        control_system::FutureMeasurements measurements(measurements_per_update,
                                                        0.0);
        measurements.update(timescale);
        while (measurements.next_measurement().value_or(
                   std::numeric_limits<double>::infinity()) <=
               initial_id.step_time().value()) {
          measurements.pop_front();
        }
        return measurements;
      };

  auto measurementsA = setup_measurements(*timescales["LabelA"]);
  auto measurementsBC = setup_measurements(*timescales["LabelBLabelC"]);

  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<Metavariables>;
  using component = Component<Metavariables>;
  MockRuntimeSystem runner{
      {}, {std::move(timescales), std::move(functions_of_time)}};
  ActionTesting::emplace_array_component_and_initialize<component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, element_id,
      {std::move(domain), std::move(measurementsA), std::move(measurementsBC),
       std::move(stepper), initial_id, TimeStepId{},
       (initial_step_end - initial_time) * slab.duration(), TimeDelta{},
       AdaptiveSteppingDiagnostics{}, std::move(history)});

  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
  const bool ready = ActionTesting::next_action_if_ready<component>(
      make_not_null(&runner), element_id);
  if (ready and expected_step_end.has_value()) {
    const TimeDelta expected_step =
        (*expected_step_end - initial_time) * slab.duration();
    CHECK(ActionTesting::get_databox_tag<component, Tags::TimeStep>(
              runner, element_id) == expected_step);
  } else {
    CHECK(not ready);
    CHECK(not expected_step_end.has_value());
  }
}

void test(const std::string& test_label, const Rational initial_time,
          const Rational initial_step_end,
          const std::optional<Rational>& expected_step_end,
          const std::vector<std::pair<double, double>>& measurement_updatesA,
          const std::vector<std::pair<double, double>>& fot_updatesA,
          const std::vector<std::pair<double, double>>& measurement_updatesBC,
          const std::vector<std::pair<double, double>>& fot_updatesB,
          const std::vector<std::pair<double, double>>& fot_updatesC,
          const std::unique_ptr<LtsTimeStepper> stepper =
              std::make_unique<TimeSteppers::AdamsMoultonPc>(3),
          TimeSteppers::History<Var::type> history = {}) {
  test_impl(test_label + ", in region", true, initial_time, initial_step_end,
            expected_step_end, measurement_updatesA, fot_updatesA,
            measurement_updatesBC, fot_updatesB, fot_updatesC,
            serialize_and_deserialize(stepper), history);
  test_impl(test_label + ", not in region", false, initial_time,
            initial_step_end, initial_step_end, measurement_updatesA,
            fot_updatesA, measurement_updatesBC, fot_updatesB, fot_updatesC,
            serialize_and_deserialize(stepper), history);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Actions.LimitTimeStepLts",
                  "[Unit][ControlSystem]") {
  register_classes_with_charm<
      TimeSteppers::AdamsBashforth, TimeSteppers::AdamsMoultonPc,
      domain::FunctionsOfTime::PiecewisePolynomial<0>,
      domain::CoordinateMap<Frame::BlockLogical, Frame::Inertial,
                            domain::CoordinateMaps::Identity<1>>,
      domain::CoordinateMap<Frame::BlockLogical, Frame::Grid,
                            domain::CoordinateMaps::Identity<1>>,
      domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                            domain::CoordinateMaps::Identity<1>>,
      domain::CoordinateMap<Frame::Grid, Frame::Distorted,
                            domain::CoordinateMaps::Identity<1>>,
      domain::CoordinateMap<Frame::Distorted, Frame::Inertial,
                            domain::CoordinateMaps::Identity<1>>>();
  const double infinity = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::signaling_NaN();
  const double arbitrary = -1234.0;

  // FIXME goal interval not relevant?
  // For each active control system below, the interval from the
  // measurement triggering the update to the FoT expiration time is
  // indicated as the goal interval.  A step must occur in this
  // interval to avoid a deadlock.

  // FIXME should adjust values to stay in slab
  // clang-format off
  test("No control systems", {1, 2}, {3, 4}, {{3, 4}},
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}});
  test("Step much shorter than limit", {1, 2}, {3, 4}, {{3, 4}},
       {{0.0, 1.0}, {5.0, nan}},
       {{0.0, arbitrary}, {5.0, nan}},  // goal range [20, 50]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}});
  test("Limited by expiration", {1, 2}, {3, 4}, {{5, 8}},
       {{0.0, 0.3}, {1.0, nan}},
       {{0.0, arbitrary}, {0.7, nan}},  // goal range [2, 4]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}});
  // test("Adjusted to keep steps even", 1.0, 8.0, 6.0,
  //      {{0.0, 5.0}, {20.0, nan}},
  //      {{0.0, arbitrary}, {11.0, nan}},  // goal range [10, 11]
  //      {{0.0, infinity}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}});
  // test("Past update but not to expiration", 1.0, 8.0, 8.0,
  //      {{0.0, 3.0}, {20.0, nan}},
  //      {{0.0, arbitrary}, {11.0, nan}},  // goal range [6, 11]
  //      {{0.0, infinity}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}});
  test("Limited by expiration, 2 systems", {1, 2}, {3, 4}, {{5, 8}},
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, 0.3}, {1.0, nan}},
       {{0.0, arbitrary}, {0.9, nan}},  // goal range [2, 5]
       {{0.0, arbitrary}, {0.7, nan}});  // goal range [2, 6]
  test("Limited by expiration, 2 measurements", {1, 2}, {3, 4}, {{5, 8}},
       {{0.0, 0.3}, {1.0, nan}},
       {{0.0, arbitrary}, {0.7, nan}},  // goal range [4, 6]
       {{0.0, 0.4}, {1.0, nan}},
       {{0.0, arbitrary}, {0.9, nan}},  // goal range [2, 5]
       {{0.0, arbitrary}, {infinity, nan}});
  // test("Adjusted to keep steps even, 2 systems", 1.0, 11.0, 7.0,
  //      {{0.0, 6.0}, {20.0, nan}},
  //      {{0.0, arbitrary}, {13.0, nan}},  // goal range [12, 13]
  //      {{0.0, 2.0}, {20.0, nan}},
  //      {{0.0, arbitrary}, {15.0, nan}},  // goal range [4, 15]
  //      {{0.0, arbitrary}, {infinity, nan}});
  // test("Adjusted to keep steps even, limited by update", 1.0, 11.0, 8.0,
  //      {{0.0, 6.0}, {20.0, nan}},
  //      {{0.0, arbitrary}, {13.0, nan}},  // goal range [12, 13]
  //      {{0.0, 4.0}, {20.0, nan}},
  //      {{0.0, arbitrary}, {15.0, nan}},  // goal range [8, 15]
  //      {{0.0, arbitrary}, {infinity, nan}});
  test("Would be limited by expiration as of now", {1, 2}, {3, 4}, {{3, 4}},
       {{0.0, 0.3}, {1.0, nan}},
       {{0.0, arbitrary}, {0.55, arbitrary}, {2.0, nan}},  // goal range [4, 20]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}});
  test("Insufficient timescale data", {1, 2}, {3, 4}, std::nullopt,
       {{0.0, 1.0}, {0.5, nan}},
       {{0.0, arbitrary}, {5.0, nan}},  // goal range [?, ?]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}});
  // test("Barely sufficient timescale data", 1.0, 5.0, 5.0,
  //      {{0.0, 10.0}, {10.0, nan}},
  //      {{0.0, arbitrary}, {50.0, nan}},  // goal range [20, 50]
  //      {{0.0, infinity}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}});
  test("Insufficient FoT data", {1, 2}, {3, 4}, std::nullopt,
       {{0.0, 1.0}, {5.0, nan}},
       {{0.0, arbitrary}, {1.5, nan}},  // goal range [20, ?]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}});
  // test("Barely insufficient FoT data", 1.0, 5.0, std::nullopt,
  //      {{0.0, 10.0}, {50.0, nan}},
  //      {{0.0, arbitrary}, {20.0, nan}},  // goal range [20, ?]
  //      {{0.0, infinity}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}},
  //      {{0.0, arbitrary}, {infinity, nan}});

  test("Does nothing with Adams-Bashforth", {1, 2}, {3, 4}, {{3, 4}},
       {{0.0, 0.3}, {0.5, nan}},
       {{0.0, arbitrary}, {0.7, nan}},  // goal range [2, 3]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       std::make_unique<TimeSteppers::AdamsBashforth>(4));
  test("Doesn't need data with Adams-Bashforth", {1, 2}, {3, 4}, {{3, 4}},
       {{0.0, 10.0}, {5.0, nan}},
       {{0.0, arbitrary}, {50.0, nan}},  // goal range [?, ?]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       std::make_unique<TimeSteppers::AdamsBashforth>(4));

  TimeSteppers::History<Var::type> history{};
  {
    const TimeStepId now(true, 0, Time(Slab(0.0, 1.0), {1, 2}));
    const TimeStepId starting(true, -1, Time(Slab(0.0, 1.0), {3, 4}));
    history.insert(starting, {}, {});
    history.insert(now, {}, {});
    REQUIRE(
        not TimeSteppers::AdamsMoultonPc(4).can_change_step_size(now, history));
  }
  test("Step can't change but is OK", {1, 2}, {3, 4}, {{3, 4}},
       {{0.0, 10.0}, {50.0, nan}},
       {{0.0, arbitrary}, {50.0, nan}},  // goal range [20, 50]
       {{0.0, infinity}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       {{0.0, arbitrary}, {infinity, nan}},
       std::make_unique<TimeSteppers::AdamsMoultonPc>(4), history);
  CHECK_THROWS_WITH(
      test("Step can't change but is not OK", {1, 2}, {3, 4}, {{0}},
           {{0.0, 0.3}, {5.0, nan}},
           {{0.0, arbitrary}, {0.7, nan}},  // goal range [2, 4]
           {{0.0, infinity}, {infinity, nan}},
           {{0.0, arbitrary}, {infinity, nan}},
           {{0.0, arbitrary}, {infinity, nan}},
           std::make_unique<TimeSteppers::AdamsMoultonPc>(4), history),
      Catch::Matchers::ContainsSubstring(
          "Step must be decreased to avoid control-system deadlock, but "
          "time-stepper requires a fixed step size."));
  // clang-format on
}
