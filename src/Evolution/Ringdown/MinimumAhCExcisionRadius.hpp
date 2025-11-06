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
 * \brief This function finds a safe ringdown excision radius for starting
 * the ringdown of a common horizon from a binary inspiral.

 * \details It does this by constructing inspiral-grid-frame AhA/AhB excision
 * strahlkorpers from the radii and centers passed to this function. It then
 * maps those excisions to the inspiral-inertial-frame using the inspiral domain
 * and functions of time from the inspiral volume data and subfile supplied. We
 * then construct a test ringdown domain that has all the corrected functions of
 * time and an initial guess for the inner radius. The excision shape
 * coefficients are expected to be in an h5 file specified by
 * `path_to_AhC_distorted_h5` and `AhC_distorted_subfile_names`. The expansion
 * and rotation maps should be the maps used in the inspiral specified by
 * `exp_func_and_2_derivs`, `exp_outer_bdry_func_and_2_derivs`, and
 * `rot_func_and_2_derivs`. The translation map is not from the inspiral, it
 * should be the corrected function and two derivatives returned by
 * ComputeAhCCoefsInRingdownDistortedFrame.py specified through
 * `trans_func_and_2_derivs`. The inspiral-inertial-frame excision A/B are then
 * transformed to the ringdown-grid-frame using the test ringdown domain and
 * ringdown functions of time described above. The inner radius is iterated
 * upon, using 2 main loops, the outer loop that changes the L_max for the
 * excisions A/B being transformed from the inspiral and the inner loop that
 * changes the excision radius. The outer loop converges when multiple L_max
 * values for excisions A/B fit inside the proposed rindown domain and the
 * difference between the excision radius used in the previous outer loop
 * iteration and current outer loop iteration are within a tolerance set by
 * 1e-3 / q where q is the mass ratio. The inner loop converges when the
 * difference between the excision radius used in the previous inner loop
 * iteration and current inner loop iteration are within a tolerance set by
 * 1e-3 / q.
 */
namespace evolution::Ringdown {

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
