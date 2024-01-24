// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <memory>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Time/Actions/TakeStep.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/IsUsingTimeSteppingErrorControl.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/Rk3HesthavenSsp.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace PUP {
struct er;
}  // namespace PUP

namespace {
struct Var : db::SimpleTag {
  using type = double;
};

struct System {
  using variables_tag = Var;
};

template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using const_global_cache_tags = tmpl::list<>;
  using simple_tags =
      tmpl::list<Tags::ConcreteTimeStepper<TimeStepper>, Tags::TimeStep,
                 ::Tags::IsUsingTimeSteppingErrorControl, Var,
                 Tags::HistoryEvolvedVariables<Var>>;
  using compute_tags = time_stepper_ref_tags<TimeStepper>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 simple_tags, compute_tags>>>,
      Parallel::PhaseActions<Parallel::Phase::Testing,
                             tmpl::list<Actions::TakeStep<System, false>>>>;
};

struct Metavariables {
  using component_list = tmpl::list<Component<Metavariables>>;

  void pup(PUP::er& /*p*/) {}
};

SPECTRE_TEST_CASE("Unit.Time.Actions.TakeStep", "[Unit][Time][Actions]") {
  register_classes_with_charm<TimeSteppers::Rk3HesthavenSsp>();

  using metavariables = Metavariables;
  using component = Component<metavariables>;

  const Slab slab(1., 3.);
  const TimeDelta time_step = slab.duration() / 2;

  Var::type vars = 1.0;
  Tags::HistoryEvolvedVariables<Var>::type history{3};
  history.insert(TimeStepId(true, 0, slab.start()), vars, 3.0);

  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavariables>;
  MockRuntimeSystem runner{{}};
  ActionTesting::emplace_component_and_initialize<component>(
      &runner, 0,
      {std::make_unique<TimeSteppers::Rk3HesthavenSsp>(), time_step, false,
       vars, std::move(history)});

  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  auto box_for_function = serialize_and_deserialize(
      ActionTesting::get_databox<component>(make_not_null(&runner), 0));

  runner.next_action<component>(0);
  take_step<System, false>(make_not_null(&box_for_function));

  const auto& action_box = ActionTesting::get_databox<component>(runner, 0);
  CHECK(db::get<Var>(box_for_function) == db::get<Var>(action_box));
}
}  // namespace
