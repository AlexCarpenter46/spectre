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
 * \brief This function is used to transition from inspiral to ringdown. It
 * reads inertial frame common horizon Strahlkorper coefs from a file and
 * returns the Strahlkorper's ringdown grid frame coefs and the inertial frame
 * geometric centers at multiple times which will be used to initialize the
 * shape and translation function of time for the ringdown.
 *
 * \details Reads common horizon Strahlkorpers (assumed to be in the inspiral
 * inertial frame) from a file, then transforms them to a temporary ringdown
 * domain defined by the expansion, rotation, and translation maps from the
 * inspiral specified by `exp_func_and_2_derivs`,
 * `exp_outer_bdry_func_and_2_derivs`, `rot_func_and_2_derivs`, and
 * `trans_func_and_2_derivs`. The expansion and rotation functions of time
 * from the inspiral are the same as the ringdown frame's expansion and rotation
 * maps at the given `match_time`, but later settle to constant values by the
 * given `settling_timescale`. The translation function of time supplied from
 * the inspiral is not the translation map we'll use in the final ringdown, it
 * is only used to correctly map the common horizon's geometric center so that
 * we can make the corrected translation map for the final ringdown. A shape map
 * is not specified, because we do not yet know the shape coefficients of the
 * common horizon for the ringdown. We get the shape coefficients by
 * transforming the common horizon to the temporary ringdown grid frame which is
 * almost the temporary ringdown distorted frame except for the uncorrected
 * translation map. The translation map is corrected by transforming the common
 * horizon from the temporary ringdown grid frame to the temporary ringdown
 * inertial frame and saving the geometric center points at multiple times. We
 * then take those geometric center points and build the corrected translation
 * function of time. This is done because the center of the Strahlkorper in this
 * temporary ringdown grid frame is NOT the origin (this is the case whether or
 * not you say Recenter=true when transforming the Strahlkorper), but the
 * distortion map does its distortion about the origin.
 * Possible ways to account for this:
 *  1) Put a translation map before the distortion map.
 *  2) Change the center of the distortion map.  (then we need to live with
 *     this during the ringdown).
 *  3) Correct the current translation map so that the excision boundary
 *     maps to the correct place.
 * We choose 3). Section 6 of https://arxiv.org/pdf/1211.6079 explains in more
 * detail of how we initialize the shape/translation map, but idea is that in
 * Eq. (104), the horizon can be written as
 * x^ibar_AH = x^ibar_AHc + Sum(Slm Ylm) n^ibar(theta,phi) where x^ibar_AHc is
 * the center of the Strahlkorper in the temporary ringdown grid frame, which is
 * time-dependent and not zero, n^ibar is the direction unit vector in
 * the (theta,phi) direction relative to x^ibar_AHc, and ibar is the index
 * corresponding to the temporary ringdown grid frame.
 * Now x^i = T0^i + M^i_ibar x^ibar where T0^i is the inspiral translation map,
 * and M^i_ibar is scaling+rotation.
 * Thus
 * x^i_AH = T0^i + M^i_ibar x^ibar_AH
 *        = T0^i + M^i_ibar x^ibar_AHc + M^i_ibar Sum(Slm Ylm) n^ibar
 * Therefore if you define a new translation map T^i as in Eq. (107) (corrected
 * translation map that will be used in the final ringdown)
 * T^i = T0^i + M^i_ibar x^ibar_AHc (that is, you define T^i to be the
 * same as x^i_AHc), then you can rewrite the relationship as
 * x^i_AH = T^i + M^i_ibar Sum(Slm Ylm) n^ibar
 * Therefore we use a new map x^i = T^i + M^i_ibar x^itilde where itilde refers
 * to the final ringdown distorted frame, x^itilde is a new coordinate where
 * x^itilde_AHc = 0 and the coefficients of the AH in the x^itilde frame can be
 * used unchanged (except for a minus sign) in the distortion map that connects
 * x^igrid and x^itilde.
 * \note Only temporary ringdown grid frame common horizon coefs and
 * temporary ringdown inertial frame geometric center points within
 * `requested_number_of_times_from_end` times from the final time are returned.
 */
std::pair<std::vector<DataVector>, std::vector<std::array<double, 3>>>
strahlkorper_coefs_and_centers(
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
