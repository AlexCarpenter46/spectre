// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <memory>
#include <optional>

#include "ControlSystem/CombinedName.hpp"
#include "ControlSystem/Metafunctions.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/Callback.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Time/ChooseLtsStepSize.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Time.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct TimeStep;
struct TimeStepId;
template <typename StepperInterface>
struct TimeStepper;
}  // namespace Tags
class TimeStepper;
namespace control_system::Tags {
template <typename ControlSystems>
struct FutureMeasurements;
struct MeasurementTimescales;
}  // namespace control_system::Tags
namespace domain::Tags {
template <size_t VolumeDim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
namespace tuples {
template <typename... Tags>
class TaggedTuple;
}  // namespace tuples
/// \endcond

namespace control_system::Actions {
// FIXME
/// \ingroup ControlSystemsGroup
/// \brief Limit the step size in a GTS evolution to prevent deadlocks from
/// control system measurements.
///
/// \details Most time steppers require evaluations of the coordinates
/// at several times during the step before they can produce dense
/// output.  If any of those evaluations require a function-of-time
/// update depending on a measurement within the step, the evolution
/// will deadlock.  This action reduces the step size if necessary to
/// prevent that from happening.
///
/// Specifically:
/// 1. The chosen step will never be longer than the unmodified step,
///    and will be short enough to avoid relevant function-of-time
///    expirations.
/// 2. Given the previous, the step will cover as many control-system
///    updates as possible.
/// 3. If the next step is likely to be limited by this action, adjust
///    the length of the current step so that this step and the next
///    step will be as close as possible to the same size.
template <typename ControlSystems>
struct LimitTimeStepLts {
 private:
  using control_system_groups =
      tmpl::transform<metafunctions::measurements_t<ControlSystems>,
                      metafunctions::control_systems_with_measurement<
                          tmpl::pin<ControlSystems>, tmpl::_1>>;

  template <typename Group>
  struct GroupExpiration {
    using type = double;
  };

 public:
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Metavariables::volume_dim>& element_id,
      ActionList /*meta*/, const ParallelComponent* const /*meta*/) {
    static_assert(Metavariables::local_time_stepping,
                  "The control system LimitTimeStepLts action is only for "
                  "local time stepping.");

    // FIXME this should use different regions for different control systems
    const auto& block =
        db::get<domain::Tags::Domain<Metavariables::volume_dim>>(box)
            .blocks()[element_id.block_id()];
    if (not block.has_distorted_frame()) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const auto& time_step_id = db::get<::Tags::TimeStepId>(box);
    if (time_step_id.substep() != 0) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const auto& time_stepper = db::get<::Tags::TimeStepper<TimeStepper>>(box);
    if (time_stepper.number_of_substeps() == 1) {
      // If there are no substeps, there is no reason to limit the
      // step size so substeps can be evaluated.  Single-substep FSAL
      // methods could be problematic, but we don't have any of those.
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const auto& this_proxy =
        ::Parallel::get_parallel_component<ParallelComponent>(
            cache)[element_id];

    // Minimum expiration time for any FoT in the measurement group.
    tmpl::wrap<tmpl::transform<control_system_groups,
                               tmpl::bind<GroupExpiration, tmpl::_1>>,
               tuples::TaggedTuple>
        group_expiration_times{};

    bool ready = true;
    // Calculate group_expiration_times
    tmpl::for_each<control_system_groups>([&](auto group_v) {
      if (not ready) {
        return;
      }
      using group = tmpl::type_from<decltype(group_v)>;

      auto& future_measurements =
          db::get_mutable_reference<Tags::FutureMeasurements<group>>(
              make_not_null(&box));

      std::optional<double> group_update = future_measurements.next_update();
      if (not group_update.has_value()) {
        Parallel::mutable_cache_item_is_ready<
            control_system::Tags::MeasurementTimescales>(
            cache,
            Parallel::make_array_component_id<ParallelComponent>(element_id),
            [&](const auto& measurement_timescales) {
              const auto& group_timescale =
                  *measurement_timescales.at(combined_name<group>());
              future_measurements.update(group_timescale);
              group_update = future_measurements.next_update();
              ready = group_update.has_value();
              return ready ? std::unique_ptr<Parallel::Callback>{}
                           : std::unique_ptr<Parallel::Callback>(
                                 new Parallel::PerformAlgorithmCallback(
                                     this_proxy));
            });
        if (not ready) {
          return;
        }
      }

      auto& group_expiration =
          get<GroupExpiration<group>>(group_expiration_times);
      group_expiration = std::numeric_limits<double>::infinity();

      if (*group_update == std::numeric_limits<double>::infinity()) {
        // Control measurement is not active.
        return;
      }

      // Calculate group_expiration
      Parallel::mutable_cache_item_is_ready<domain::Tags::FunctionsOfTime>(
          cache,
          Parallel::make_array_component_id<ParallelComponent>(element_id),
          [&](const auto& functions_of_time) {
            tmpl::for_each<group>([&](auto system) {
              using System = tmpl::type_from<decltype(system)>;
              if (not ready) {
                return;
              }
              const auto& fot = *functions_of_time.at(System::name());
              ready = fot.time_bounds()[1] > *group_update;
              if (ready) {
                group_expiration = std::min(
                    group_expiration, fot.expiration_after(*group_update));
              }
            });
            return ready ? std::unique_ptr<Parallel::Callback>{}
                         : std::unique_ptr<Parallel::Callback>(
                               new Parallel::PerformAlgorithmCallback(
                                   this_proxy));
          });
    });
    if (not ready) {
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    // const double orig_step_start = time_step_id.step_time().value();
    const double orig_step_end =
        (time_step_id.step_time() + db::get<::Tags::TimeStep>(box)).value();

    // Smallest of the current step end and the FoT expirations.  We
    // can't step any farther than this.
    const double latest_valid_step =
        tmpl::as_pack<control_system_groups>([&](auto... groups) {
          return std::min(
              {orig_step_end,
               get<GroupExpiration<tmpl::type_from<decltype(groups)>>>(
                   group_expiration_times)...});
        });

    // if (not time_stepper.can_change_step_size(
    //         time_step_id, db::get<::Tags::HistoryEvolvedVariables<>>(box))) {
    //   if (orig_step_end > latest_valid_step) {
    //     ERROR(
    //         "Step must be decreased to avoid control-system deadlock, but "
    //         "time-stepper requires a fixed step size.");
    //   }
    //   return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    // }

#if 0
    // The last update that we can perform on the next step.  Don't
    // shrink the step past this time since that will force another
    // step to take the measurement.
    double last_update_time = orig_step_start;
    // Step time that produces a balanced step with the following
    // step, ignoring the restrictions on the current step.
    double preferred_step_time = orig_step_end;

    tmpl::for_each<control_system_groups>([&](auto group_v) {
      using group = tmpl::type_from<decltype(group_v)>;

      // This was used above, so it is not nullopt.
      const double group_update =
          db::get<Tags::FutureMeasurements<group>>(box).next_update().value();
      if (group_update <= latest_valid_step) {
        // We've satisfied this measurement.
        last_update_time = std::max(last_update_time, group_update);
      } else {
        // We can't make it far enough to do the final measurement.
        // Try to avoid a small step by choosing two equal-sized steps
        // to the expiration time.
        const double equal_step_time =
            0.5 * (orig_step_start +
                   get<GroupExpiration<group>>(group_expiration_times));
        preferred_step_time = std::min(preferred_step_time, equal_step_time);
      }
    });
#endif

    // FIXME what should we use here?
    // GTS does:
    // const double new_step_end =
    //     std::clamp(preferred_step_time, last_update_time, latest_valid_step);
    const TimeDelta new_step_size = choose_lts_step_size(
        time_step_id.step_time(),
        // FIXME refactor so step size is variable
        latest_valid_step - time_step_id.step_time().value());
    if (new_step_size == db::get<::Tags::TimeStep>(box)) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    if (not time_stepper.can_change_step_size(
            time_step_id, db::get<::Tags::HistoryEvolvedVariables<>>(box))) {
      ERROR(
          "Step must be decreased to avoid control-system deadlock, but "
          "time-stepper requires a fixed step size.");
    }

    db::mutate<::Tags::Next<::Tags::TimeStepId>, ::Tags::TimeStep>(
        [&new_step_size, &time_step_id, &time_stepper](
            const gsl::not_null<TimeStepId*> local_next_time_id,
            const gsl::not_null<TimeDelta*> local_time_step) {
          *local_time_step = new_step_size;
          *local_next_time_id =
              time_stepper.next_time_id(time_step_id, new_step_size);
        },
        make_not_null(&box));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace control_system::Actions
