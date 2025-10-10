// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"

/*!
 * \brief Functionality for evolving a ringdown following a compact-binary
 * merger.
 */
namespace evolution::Ringdown {
/*!
 * \brief Transforms `Strahlkorper` coefs to ringdown distorted frame and tracks
 * the center of mass.
 *
 * \details Reads Strahlkorper coefficients (assumed to be in the inertial
 * frame) from a file, then transforms them into the
 * ringdown distorted frame defined by the expansion, rotation, and translation
 * maps specified by `exp_func_and_2_derivs`,
 * `exp_outer_bdry_func_and_2_derivs`, `rot_func_and_2_derivs`, and
 * `trans_func_and_2_derivs`. The expansion and rotation functions of time
 * correspond to the ringdown frame's expansion and rotation maps at the time
 * given by `match_time`, and by `settling_timescale`, the timescale for the
 * maps to settle to constant values. The translation function of time supplied
 * does not correspond to the ringdown frame's translation function of time. The
 * ringdown's translation function of time needs to be built by tracking the
 * position of the center of mass at multiple times. Only Strahlkorpers and
 * center of mass points within `requested_number_of_times_from_end` times from
 * the final time are returned. This function is used to transition from
 * inspiral to ringdown; in this case, the inertial-frame Strahlkorper is the
 * common apparent horizon from a binary-black-hole inspiral; the
 * ringdown-distorted-frame coefficients are used to initialize the shape map
 * for the ringdown domain. The center of mass points are used to initialize the
 * translation map for the ringdown domain.
 */
std::pair<std::vector<DataVector>, std::vector<std::array<double, 3>>>
strahlkorper_coefs_in_ringdown_distorted_frame(
    const std::string& path_to_volume_data,
    const std::string& volume_subfile_name,
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name,
    size_t requested_number_of_times_from_end, double match_time,
    double settling_timescale,
    const std::optional<std::array<double, 3>>& exp_func_and_2_derivs =
        std::nullopt,
    const std::optional<std::array<double, 3>>&
        exp_outer_bdry_func_and_2_derivs = std::nullopt,
    const std::optional<std::vector<std::array<double, 4>>>&
        rot_func_and_2_derivs = std::nullopt,
    const std::optional<std::array<std::array<double, 3>, 3>>&
        trans_func_and_2_derivs = std::nullopt);
}  // namespace evolution::Ringdown
