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
 * \brief This function finds the minimum radius excision factor needed
 * for continuing evolutions into ringdown.
 */
namespace evolution::Ringdown {

// What I need so far: AhC horizon file path, subfile name, shape coefficients
// from computeAhCCeofsinRingdownDistortedFrame,
// requested_number_of_times_from_end, match_time, settling_timescale, fots,
// radius_excision_A, radius_excision_B, excision_center_A, excision_center_B

// I think that's everything?

double minimum_ahc_excision_radius(
    const std::string& path_to_volume_data,
    const std::string& volume_subfile_name,
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name,
    const std::string& path_to_AhC_distorted_h5,
    const std::vector<std::string>& AhC_distorted_subfile_names,
    size_t requested_number_of_times_from_end, double match_time,
    double settling_timescale, double excision_A_radius,
    double excision_B_radius, std::array<double, 3> excision_A_center,
    std::array<double, 3> excision_B_center, size_t excision_l_max,
    const std::optional<std::array<double, 3>>& exp_func_and_2_derivs =
        std::nullopt,
    const std::optional<std::array<double, 3>>&
        exp_outer_bdry_func_and_2_derivs = std::nullopt,
    const std::optional<std::vector<std::array<double, 4>>>&
        rot_func_and_2_derivs = std::nullopt,
    const std::optional<std::array<std::array<double, 3>, 3>>&
        trans_func_and_2_derivs = std::nullopt);
}  // namespace evolution::Ringdown
