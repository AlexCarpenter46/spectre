// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"

/*!
 * \brief This function finds the minimum safe ahc excision radius needed
 * for starting the ringdown of a common horizon. It does this by taking the
 * inspiral AhA/AhB excision strahlkorpers and transforms them to a ringdown
 * domain containing all of the correct time dependent maps for the ringdown.
 * There are 2 main loops, the outer loop that changes the excision radius and
 * the inner loop that changes the l_max for the excisions A/B being transformed
 * from the inspiral. The inner loop converges when multiple l_max values for
 * excisions A/B fit inside the proposed rindown domain. The outer loop
 * converges when the difference between the excision radius used in the
 * previous iteration and current iteration are within a tolerance set by 1e-3 /
 * q where q is the mass ratio.
 */
namespace evolution::Ringdown {

std::pair<double, double> minimum_ahc_excision_radius(
    const std::string& path_to_volume_data,
    const std::string& volume_subfile_name,
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name,
    const std::string& path_to_AhC_distorted_h5,
    const std::vector<std::string>& AhC_distorted_subfile_names,
    size_t requested_number_of_times_from_end, double match_time,
    double settling_timescale, double excision_A_radius,
    double excision_B_radius, std::array<double, 3> excision_A_center,
    std::array<double, 3> excision_B_center,
    const std::optional<std::array<double, 3>>& exp_func_and_2_derivs =
        std::nullopt,
    const std::optional<std::array<double, 3>>&
        exp_outer_bdry_func_and_2_derivs = std::nullopt,
    const std::optional<std::vector<std::array<double, 4>>>&
        rot_func_and_2_derivs = std::nullopt,
    const std::optional<std::array<std::array<double, 3>, 3>>&
        trans_func_and_2_derivs = std::nullopt);
}  // namespace evolution::Ringdown
