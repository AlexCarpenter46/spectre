// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Time/Actions/ChangeStepSize.hpp"
#include "Time/AdaptiveSteppingDiagnostics.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/PreviousStepperError.hpp"
#include "Time/Tags/StepperError.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
struct AllStepChoosers;
namespace Tags {
struct AdaptiveSteppingDiagnostics;
struct IsUsingTimeSteppingErrorControl;
struct StepperErrorUpdated;
struct TimeStep;
template <typename StepperInterface>
struct TimeStepper;
}  // namespace Tags
/// \endcond

namespace take_step_detail {
template <typename System, typename VariablesTag, typename DbTags>
void update_one_variables(const gsl::not_null<db::DataBox<DbTags>*> box) {
  using history_tag = Tags::HistoryEvolvedVariables<VariablesTag>;
  bool is_using_error_control = false;
  if constexpr (db::tag_is_retrievable_v<Tags::IsUsingTimeSteppingErrorControl,
                                         db::DataBox<DbTags>>) {
    is_using_error_control =
        db::get<Tags::IsUsingTimeSteppingErrorControl>(*box);
  }
  if (is_using_error_control) {
    using error_tag = ::Tags::StepperError<VariablesTag>;
    using previous_error_tag = ::Tags::PreviousStepperError<VariablesTag>;
    if constexpr (tmpl::list_contains_v<DbTags, error_tag>) {
      db::mutate<Tags::StepperErrorUpdated, VariablesTag, error_tag,
                 previous_error_tag, history_tag>(
          [](const gsl::not_null<bool*> stepper_error_updated,
             const gsl::not_null<typename VariablesTag::type*> vars,
             const gsl::not_null<typename error_tag::type*> error,
             const gsl::not_null<typename previous_error_tag::type*>
                 previous_error,
             const gsl::not_null<typename history_tag::type*> history,
             const ::TimeDelta& time_step, const TimeStepper& time_stepper) {
            using std::swap;
            set_number_of_grid_points(previous_error, *vars);
            swap(*error, *previous_error);
            *stepper_error_updated = time_stepper.update_u(
                vars, make_not_null(&*error), history, time_step);
            if (not *stepper_error_updated) {
              swap(*error, *previous_error);
            }
          },
          box, db::get<Tags::TimeStep>(*box),
          db::get<Tags::TimeStepper<TimeStepper>>(*box));
    } else {
      ERROR(
          "Cannot update the stepper error measure -- "
          "`::Tags::StepperError<VariablesTag>` is not present in the box.");
    }
  } else {
    db::mutate<VariablesTag, history_tag>(
        [](const gsl::not_null<typename VariablesTag::type*> vars,
           const gsl::not_null<typename history_tag::type*> history,
           const ::TimeDelta& time_step, const TimeStepper& time_stepper) {
          time_stepper.update_u(vars, history, time_step);
        },
        box, db::get<Tags::TimeStep>(*box),
        db::get<Tags::TimeStepper<TimeStepper>>(*box));
  }
}

template <typename System, typename DbTags>
void update_variables(const gsl::not_null<db::DataBox<DbTags>*> box) {
  if constexpr (tt::is_a_v<tmpl::list, typename System::variables_tag>) {
    // The system has multiple evolved variables, probably because
    // there is a mixture of real and complex values or similar.  Step
    // all of them.
    tmpl::for_each<typename System::variables_tag>([&](auto tag) {
      update_one_variables<System, tmpl::type_from<decltype(tag)>>(box);
    });
  } else {
    update_one_variables<System, typename System::variables_tag>(box);
  }
}
}  // namespace take_step_detail

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
      take_step_detail::update_variables<System>(box);
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
    take_step_detail::update_variables<System>(box);
  }
}
