// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"

namespace evolution::Ringdown {

/*!
 * \brief This function finds a safe ringdown excision radius for starting
 * the ringdown of a common horizon from a binary inspiral. It does this by
 * finding a radius that will enclose every point on strahlkorpers that
 * represent excisions A/B from the inspiral-inertial-frame at the match time.
 * It does this by taking excisions A/B from the inspiral-inertial-frame and
 * transforming them to the ringdown-grid-frame. It iterates over all the points
 * from excisions A/B in the ringdown-grid-frame until it finds the point that's
 * the farthest from the geometric center of the ringdown excision. This is now
 * the minimum radius we can choose to enclose the excisions A/B from the
 * inspiral. This is not the ideal radius however, we then choose a radius
 * that's 3/4 of the way between the radius of the common horizon at the match
 * time and the minimum radius that encloses excisions A/B from the inspiral.

 * \details It does this by constructing inspiral-grid-frame AhA/AhB excision
 * strahlkorpers from the radii and centers passed to this function. It then
 * maps those excisions to the inspiral-inertial-frame using the inspiral domain
 * and functions of time from the inspiral volume data and subfile supplied. We
 * then construct a test ringdown domain that has all the corrected functions of
 * time and an initial guess for the inner radius.
 * \param path_to_AhC_distorted_h5 Path to h5 file containing ringdown shape
 * coefficients computed using ComputeAhCCoefsInRingdownDistortedFrame.py
 * \param AhC_distorted_subfile_names Subfiles in the h5 file containing
 * shape coefficients
 * \param exp_func_and_2_derivs Expansion FoT from the inspiral
 * \param exp_outer_bdry_func_and_2_derivs Outer boundary expansion FoT from the
 * inspiral
 * \param rot_func_and_2_derivs Rotation FoT from the inspiral
 * \param trans_func_and_2_derivs Translation FoT not from the inspiral, but
 * the corrected translation FoT from ComputeAhCCoefsInRingdownDistortedFrame.py
 *
 * Using these FoTs, the inspiral-inertial-frame excisions A/B are then
 * transformed to the ringdown-grid-frame. The inner radius is iterated
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
 *
 * \note This implementation does not do everything SpEC does yet, it is
 * currently missing rescaling the shape coefficients held in
 * 'path_to_AhC_distorted_h5' by excision_radius / average_ahc_radius every
 * time it constructs a test ringdown domain. This step helps the shape of the
 * excision match the shape of the apparent horizon.
 */
double minimum_ahc_excision_radius(
    const std::string& path_to_volume_data,
    const std::string& volume_subfile_name,
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name,
    const std::string& path_to_AhC_distorted_h5,
    const std::vector<std::string>& AhC_distorted_subfile_names,
    double match_time, double settling_timescale, double excision_A_radius,
    double excision_B_radius, std::array<double, 3> excision_A_center,
    std::array<double, 3> excision_B_center,
    const std::optional<std::array<double, 3>>& exp_func_and_2_derivs =
        std::nullopt,
    const std::optional<std::array<double, 3>>&
        exp_outer_bdry_func_and_2_derivs = std::nullopt,
    const std::optional<std::vector<std::array<double, 4>>>&
        rot_func_and_2_derivs = std::nullopt,
    const std::optional<std::array<std::array<double, 3>, 3>>&
        trans_func_and_2_derivs = std::nullopt,
    double match_time_tol = 1e-12);
}  // namespace evolution::Ringdown
