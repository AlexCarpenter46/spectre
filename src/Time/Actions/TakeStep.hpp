// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>

#include "Parallel/AlgorithmExecution.hpp"
#include "Time/TakeStep.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace db {
template <typename TagsList>
class DataBox;
}  // namespace db
namespace tuples {
template <class... Tags>
class TaggedTuple;
}  // namespace tuples
/// \endcond

namespace Actions {
/// \ingroup ActionsGroup
/// \ingroup TimeGroup
/// \brief Update the evolved variables, possibly adjusting the step
/// size in LTS mode.
///
/// Uses:
/// - DataBox:
///   - system::variables_tag
///   - Tags::HistoryEvolvedVariables<variables_tag>
///   - Tags::TimeStep
///   - Tags::TimeStepper<TimeStepper>
///   - Tags::IsUsingTimeSteppingErrorControl
///   - as required for step adjustment
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies:
///   - variables_tag
///   - Tags::HistoryEvolvedVariables<variables_tag>
template <typename System, bool AdjustStepSize>
struct TakeStep {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    take_step<System, AdjustStepSize>(make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Actions
