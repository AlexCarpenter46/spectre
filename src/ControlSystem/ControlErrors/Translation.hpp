// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/math/quaternion.hpp>
#include <cstddef>
#include <iostream>
#include <pup.h>

#include "ControlSystem/ControlErrors/Expansion.hpp"
#include "ControlSystem/ControlErrors/Rotation.hpp"
#include "ControlSystem/DataVectorHelpers.hpp"
#include "ControlSystem/Protocols/ControlError.hpp"
#include "ControlSystem/Tags/QueueTags.hpp"
#include "ControlSystem/Tags/SystemTags.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/Tags/ObjectCenter.hpp"
#include "Domain/FunctionsOfTime/QuaternionHelpers.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
/// \endcond

namespace control_system {
namespace ControlErrors {
/*!
 * \brief Control error in the 3D \link
 * domain::CoordinateMaps::TimeDependent::Translation Translation \endlink
 * coordinate map
 *
 * \details Computes the error in how much the system has translated by using
 * Eq. (42) from \cite Ossokine2013zga. The equation is
 *
 * \f[ \left(0, \delta\vec{T}\right) = a\mathbf{q}\left(\mathbf{x}_A -
 * \mathbf{c}_A - \mathbf{\delta q}\wedge\mathbf{c}_A - \frac{\delta
 * a}{a}\mathbf{c}_A \right)\mathbf{q}^*
 * \f]
 *
 * where object A is located on the positive x-axis in the grid frame, bold face
 * letters are quaternions, vectors are promoted to quaternions as \f$
 * \mathbf{v} = (0, \vec{v}) \f$, \f$ \mathbf{q} \f$ is the quaternion from the
 * \link domain::CoordinateMaps::TimeDependent::Rotation Rotation \endlink map,
 * \f$ a \f$ is the function \f$ a(t) \f$ from the \link
 * domain::CoordinateMaps::TimeDependent::CubicScale CubicScale \endlink map,
 * \f$ \mathbf{\delta q}\wedge\mathbf{c}_A \equiv (0, \delta\vec{q} \times
 * \vec{c}_A) \f$, \f$ \delta\vec{q} \f$ is the \link
 * control_system::ControlErrors::Rotation Rotation \endlink control error, and
 * \f$ \delta a\f$ is the \link control_system::ControlErrors::Expansion
 * Expansion \endlink control error.
 *
 * Requirements:
 * - This control error requires that there be exactly two objects in the
 *   simulation
 * - Currently both these objects must be black holes
 * - Currently this control error can only be used with the \link
 *   control_system::Systems::Translation Translation \endlink control system
 * - There must exist an expansion map and a quaternion rotation map in the
 *   coordinate map with names "Expansion" and "Rotation", respectively.
 */
struct Translation : tt::ConformsTo<protocols::ControlError> {
  static constexpr size_t expected_number_of_excisions = 2;

  using object_centers =
      domain::object_list<domain::ObjectLabel::A, domain::ObjectLabel::B>;

  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Computes the control error for translation control. This should not "
      "take any options."};

  void pup(PUP::er& /*p*/) {}

  template <typename Metavariables, typename... TupleTags>
  DataVector operator()(const ::TimescaleTuner<true>& tuner,
                        const Parallel::GlobalCache<Metavariables>& cache,
                        const double time,
                        const std::string& /*function_of_time_name*/,
                        const tuples::TaggedTuple<TupleTags...>& measurements) {
    // if (time < 50.0) {
    //     return {0.0, 0.0, 0.0};
    // } else {
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);

    using quat = boost::math::quaternion<double>;
    std::cout << "rotation function of time: "
              << functions_of_time.at("Rotation")->func(time)[0] << std::endl;

    const quat quaternion = datavector_to_quaternion(
        functions_of_time.at("Rotation")->func(time)[0]);
    const double expansion_factor =
        functions_of_time.at("Expansion")->func(time)[0][0];

    using center_A =
        control_system::QueueTags::Center<::domain::ObjectLabel::A>;
    using center_B =
        control_system::QueueTags::Center<::domain::ObjectLabel::B>;

    const tnsr::I<double, 3, Frame::Grid>& grid_position_of_A_tnsr =
        Parallel::get<domain::Tags::ObjectCenter<domain::ObjectLabel::A>>(
            cache);
    const DataVector grid_position_of_A{{grid_position_of_A_tnsr[0],
                                         grid_position_of_A_tnsr[1],
                                         grid_position_of_A_tnsr[2]}};
    const DataVector& current_position_of_A = get<center_A>(measurements);

    const tnsr::I<double, 3, Frame::Grid>& grid_position_of_B_tnsr =
        Parallel::get<domain::Tags::ObjectCenter<domain::ObjectLabel::B>>(
            cache);
    const DataVector grid_position_of_B{{grid_position_of_B_tnsr[0],
                                         grid_position_of_B_tnsr[1],
                                         grid_position_of_B_tnsr[2]}};
    const DataVector& current_position_of_B = get<center_B>(measurements);

    const DataVector& grid_position_average =
        0.5 * (grid_position_of_A + grid_position_of_B);
    const DataVector& current_position_average =
        0.5 * (current_position_of_A + current_position_of_B);

    const DataVector& grid_separation = grid_position_of_A - grid_position_of_B;
    const DataVector& current_separation =
        current_position_of_A - current_position_of_B;

    // pMdotcM = current_separation_dot_grid_separation
    // pMdotcP = current_separation_dot_grid_average
    // cMdotcP = grid_separation_dot_grid_average
    // cMdotcM = grid_separation_dot_grid_separation
    double current_separation_dot_grid_separation = 0.0;
    double current_separation_dot_grid_average = 0.0;
    double grid_separation_dot_grid_average = 0.0;
    double grid_separation_dot_grid_separation = 0.0;
    for (size_t i = 0; i < 3; i++) {
      current_separation_dot_grid_separation +=
          current_separation[i] * grid_separation[i];
      current_separation_dot_grid_average +=
          current_separation[i] * grid_position_average[i];
      grid_separation_dot_grid_average +=
          grid_separation[i] * grid_position_average[i];
      grid_separation_dot_grid_separation +=
          grid_separation[i] * grid_separation[i];
    }

    const DataVector rotation_error =
        rotation_control_error_(tuner, cache, time, "Rotation", measurements);
    // Use A because it's on the positive x-axis, however, B would work as well.
    // Just so long as we are consistent.
    const DataVector rotation_error_cross_grid_pos_A =
        cross(rotation_error, grid_position_of_A);
    const DataVector rotation_error_cross_grid_pos_average =
        cross(rotation_error, grid_position_average);

    const double expansion_error = expansion_control_error_(
        tuner, cache, time, "Expansion", measurements)[0];

    DataVector translation_control =
        make_with_value<DataVector>(grid_position_of_A, 0.0);
    for (size_t i = 0; i < 3; i++) {
      translation_control[i] =
          expansion_factor *
          (grid_separation_dot_grid_separation * current_position_average[i] -
           current_separation_dot_grid_separation * grid_position_average[i] -
           grid_separation_dot_grid_average * current_separation[i] +
           current_separation_dot_grid_average * grid_separation[i]) /
          grid_separation_dot_grid_separation;
    }
    // From eq. 42 in 1304.3067
    // const quat middle_expression = datavector_to_quaternion(
    //     current_position_of_A -
    //     (1.0 + expansion_error / expansion_factor) * grid_position_of_A -
    //     rotation_error_cross_grid_pos_A);
    const quat middle_expression_2 = datavector_to_quaternion(
        current_position_average -
        (1.0 + expansion_error / expansion_factor) * grid_position_average -
        rotation_error_cross_grid_pos_average);
    const quat middle_expression =
        datavector_to_quaternion(translation_control);
    std::cout << "Excision separation: " << current_separation << std::endl;
    std::cout << "AH separation: " << grid_separation << std::endl;
    std::cout << "Excision average: " << current_position_average << std::endl;
    std::cout << "AH average: " << grid_position_average << std::endl;
    std::cout << "Quaternion: " << quaternion << std::endl;
    std::cout << "SpEC way: " << middle_expression << std::endl;
    std::cout << "Average/control error version: " << middle_expression_2
              << std::endl;

    // Because we are converting from a quaternion to a DataVector, there will
    // be four components in the DataVector. However, translation control only
    // requires three (the latter three to be exact, because the first component
    // should be 0. We ASSERT this also.)
    const DataVector result_with_four_components =
        expansion_factor *
        quaternion_to_datavector(quaternion * middle_expression *
                                 conj(quaternion));
    ASSERT(equal_within_roundoff(result_with_four_components[0], 0.0),
           "Error in computing translation control error. First component of "
           "resulting quaternion should be 0.0, but is " +
               get_output(result_with_four_components[0]) + " instead.");

    return {result_with_four_components[1], result_with_four_components[2],
            result_with_four_components[3]};
    //   }
  }

 private:
  Rotation rotation_control_error_{};
  Expansion expansion_control_error_{};
};
}  // namespace ControlErrors
}  // namespace control_system
