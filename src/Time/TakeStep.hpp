// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <type_traits>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Time/Actions/ChangeStepSize.hpp"
#include "Time/Actions/UpdateU.hpp"
#include "Time/AdaptiveSteppingDiagnostics.hpp"
#include "Time/Tags/AdaptiveSteppingDiagnostics.hpp"
#include "Time/Time.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
namespace Parallel::Tags {
struct Metavariables;
}  // namespace Parallel::Tags
namespace Tags {
struct TimeStep;
}  // namespace Tags
/// \endcond

/// Update the evolved variables for one substep, possibly adjusting
/// the step size in LTS mode.
template <typename System, bool AdjustStepSize,
          typename StepChoosersToUse = AllStepChoosers, typename DbTags>
void take_step(const gsl::not_null<db::DataBox<DbTags>*> box) {
  if constexpr (AdjustStepSize) {
    uint64_t step_attempts = 0;
    const auto original_step = db::get<Tags::TimeStep>(*box);
    do {
      ++step_attempts;
      update_u<System>(box);
    } while (not change_step_size<StepChoosersToUse>(box));
    db::mutate<Tags::AdaptiveSteppingDiagnostics>(
        [&](const gsl::not_null<AdaptiveSteppingDiagnostics*> diags,
            const TimeDelta& new_step) {
          diags->number_of_step_rejections += step_attempts - 1;
          if (original_step != new_step) {
            ++diags->number_of_step_fraction_changes;
          }
        },
        box, db::get<Tags::TimeStep>(*box));
  } else {
    update_u<System>(box);
  }
}
