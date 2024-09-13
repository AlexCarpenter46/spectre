# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
import warnings
from pathlib import Path
from typing import Optional, Union

import click
import numpy as np
import yaml
from rich.pretty import pretty_repr

from spectre.DataStructures import ModalVector
from spectre.SphericalHarmonics import Strahlkorper, ylm_legend_and_data


def cubic(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d


def dt_cubic(x, a, b, c, d):
    return 3 * a * x**2 + 2 * b * x + c


def dt2_cubic(x, a, b, c, d):
    return 6 * a * x + 2 * b


# Cubic fit transformed coefs to get first and second time derivatives
def fit_to_a_cubic(times, coefs, match_time, zero_coefs):
    fits = []
    fit_ahc = []
    fit_dt_ahc = []
    fit_dt2_ahc = []
    for j in np.arange(0, coefs.shape[-1], 1):
        # Optionally, avoid fitting coefficients sufficiently close to zero by
        # just setting these coefficients and their time derivatives to zero.
        if zero_coefs is not None and sum(np.abs(coefs[:, j])) < zero_coefs:
            fits.append(np.zeros(4))
            fit_ahc.append(0.0)
            fit_dt_ahc.append(0.0)
            fit_dt2_ahc.append(0.0)
            continue
        # Ignore RankWarnings suggesting the fit might not be good enough;
        # for equal-mass non-spinning, sufficiently good fits for starting
        # a ringdown, even though RankWarnings were triggered
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", np.RankWarning)
            fit = np.polyfit(times, coefs[:, j], 3)
        fits.append(fit)
        fit_ahc.append(cubic(match_time, *(fit)))
        fit_dt_ahc.append(dt_cubic(match_time, *(fit)))
        fit_dt2_ahc.append(dt2_cubic(match_time, *(fit)))

    return fit_ahc, fit_dt_ahc, fit_dt2_ahc


"""Function for Computing AhC Coefficients in the Ringdown Distorted Frame
Todo: Add description and arguments list.
"""


# Todo: Change the arguments to be h5 files containing so that we don't have to
# what's done in line 111 or 112
def compute_ahc_coefs_in_ringdown_distorted_frame(
    ahc_reductions_path,
    ahc_subfile,
    fot_vol_h5_path,
    fot_vol_subfile,
    path_to_output_h5,
    output_subfile_prefix,
    number_of_steps,
    which_obs_id,
    settling_timescale,
    zero_coefs,
):
    """Computes AhC Coefficients in the Ringdown Distorted Frame from functions
    of time from AhC in inspiral distorted frame.

    Arugments:
    ahc_reductions_path: Path to reduction file where AhC Coefficients will be
    written.
    ahc_subfile: The subfile of the reductions file where AhC coefficients will
    be placed.
    path_to_output_h5: Path to file where AhC coefficients in Ringdown distorted
    frame will be written.
    output_subfile_prefix: Subfile in output_h5 for AhC coefficients in Ringdown
    distorted frame.
    number_of_steps: The number of steps from the last time in the
    simulation to look for AhC finds.
    which_obs_id: First time to look at AhC coefficients.
    settling_timescale: Timescale for settle to constant functions of time.
    zero_coefs: Coefficients to zero out

    """

    shape_coefs = {}

    # Subfile for AhC data
    if ahc_subfile.split(".")[-1] == "dat":
        ahc_subfile = ahc_subfile.split(".")[0]

    # Setting up output with Distorted Subfile
    output_subfile_ahc = output_subfile_prefix + "AhC_Ylm"
    output_subfile_dt_ahc = output_subfile_prefix + "dtAhC_Ylm"
    output_subfile_dt2_ahc = output_subfile_prefix + "dt2AhC_Ylm"

    ahc_times = []
    ahc_legend = ""
    ahc_center = []
    ahc_lmax = 0
    with spectre_h5.H5File(ahc_reductions_path, "r") as h5file:
        datfile = h5file.get_dat("/" + ahc_subfile)
        datfile_np = np.array(datfile.get_data())
        ahc_times = datfile_np[:, 0]
        ahc_legend = datfile.get_legend()
        ahc_center = [datfile_np[0][1], datfile_np[0][2], datfile_np[0][3]]
        ahc_lmax = int(datfile_np[0][4])

    # dict storing other info for starting a ringdown obtained in the process
    # of getting the ringdown-distorted-frame AhC coefficients
    info_for_ringdown = {}
    legend_for_ringdown = {}

    # This can go somewhere else :)
    # Transform AhC coefs to ringdown distorted frame and get other data
    # needed to start a ringdown, such as initial values for functions of time
    with spectre_h5.H5File(fot_vol_h5_path, "r") as h5file:
        if fot_vol_subfile.split(".")[-1] == "vol":
            fot_vol_subfile = fot_vol_subfile.split(".")[0]
        volfile = h5file.get_vol("/" + fot_vol_subfile)
        obs_ids = volfile.list_observation_ids()
        logger.info("About to deserialize functions of time")
        fot_times = list(map(volfile.get_observation_value, obs_ids))
        serialized_fots = volfile.get_functions_of_time(obs_ids[which_obs_id])
        functions_of_time = deserialize_functions_of_time(serialized_fots)
        logger.info("Deserialized functions of time")

        # The inspiral expansion map includes two functions of time:
        # an expansion map allowing the black holes to move closer together in
        # comoving coordinates, and an expansion map causing the outer boundary
        # to move slightly inwards. The ringdown only requires the outer
        # boundary expansion map, so set the other map to the identity.
        exp_func_and_2_derivs = [1.0, 0.0, 0.0]

        exp_outer_bdry_func_and_2_derivs = [
            x[0]
            for x in functions_of_time[
                "ExpansionOuterBoundary"
            ].func_and_2_derivs(fot_times[which_obs_id])
        ]
        rot_func_and_2_derivs_tuple = functions_of_time[
            "Rotation"
        ].func_and_2_derivs(fot_times[which_obs_id])
        rot_func_and_2_derivs = [
            [coef for coef in x] for x in rot_func_and_2_derivs_tuple
        ]

        match_time = fot_times[which_obs_id]

        coefs_at_different_times = np.array(
            Ringdown.strahlkorper_coefs_in_ringdown_distorted_frame(
                ahc_reductions_path,
                ahc_subfile,
                number_of_steps,
                match_time,
                settling_timescale,
                exp_func_and_2_derivs,
                exp_outer_bdry_func_and_2_derivs,
                rot_func_and_2_derivs,
            )
        )

        # Print out coefficients for insertion into BBH domain
        logger.info("Expansion: " + str(exp_func_and_2_derivs))
        logger.info(
            "ExpansionOutrBdry: " + str(exp_outer_bdry_func_and_2_derivs)
        )
        logger.info("Rotation: " + str(rot_func_and_2_derivs))
        logger.info("Match time: " + str(match_time))
        logger.info("Settling timescale: " + str(settling_timescale))
        logger.info("Lmax: " + str(ahc_lmax))

        shape_coefs["Expansion"]: exp_func_and_2_derivs
        shape_coefs["ExpansionOutrBdry"]: exp_outer_bdry_func_and_2_derivs
        shape_coefs["Rotation"]: rot_func_and_2_derivs
        shape_coefs["Match time"]: match_time
        shape_coefs["Settling timescale"]: settling_timescale
        shape_coefs["Lmax"]: ahc_lmax

        legend_fot = ["FunctionOfTime", "dtFunctionOfTime", "dt2FunctionOfTime"]
        legend_value = ["Value"]
        info_for_ringdown["Expansion"] = exp_func_and_2_derivs
        legend_for_ringdown["Expansion"] = legend_fot
        info_for_ringdown["ExpansionOuterBdry"] = (
            exp_outer_bdry_func_and_2_derivs
        )
        legend_for_ringdown["ExpansionOuterBdry"] = legend_fot
        info_for_ringdown["Rotation0"] = [
            rot_func_and_2_derivs[0][0],
            rot_func_and_2_derivs[1][0],
            rot_func_and_2_derivs[2][0],
        ]
        info_for_ringdown["Rotation1"] = [
            rot_func_and_2_derivs[0][1],
            rot_func_and_2_derivs[1][1],
            rot_func_and_2_derivs[2][1],
        ]
        info_for_ringdown["Rotation2"] = [
            rot_func_and_2_derivs[0][2],
            rot_func_and_2_derivs[1][2],
            rot_func_and_2_derivs[2][2],
        ]
        info_for_ringdown["Rotation3"] = [
            rot_func_and_2_derivs[0][3],
            rot_func_and_2_derivs[1][3],
            rot_func_and_2_derivs[2][3],
        ]
        legend_for_ringdown["Rotation0"] = legend_fot
        legend_for_ringdown["Rotation1"] = legend_fot
        legend_for_ringdown["Rotation2"] = legend_fot
        legend_for_ringdown["Rotation3"] = legend_fot
        info_for_ringdown["MatchTime"] = [match_time]
        legend_for_ringdown["MatchTime"] = legend_value
        info_for_ringdown["SettlingTimescale"] = [settling_timescale]
        legend_for_ringdown["SettlingTimescale"] = legend_value
        info_for_ringdown["Lmax"] = [ahc_lmax]
        legend_for_ringdown["Lmax"] = legend_value

    # Do not include AhCs at times greater than the match time. Errors tend
    # to grow as time increases, so fit derivatives using the match time
    # and earlier times, to get a more accurate fit.
    ahc_times_for_fit_list = []
    coefs_at_different_times_for_fit_list = []
    for i, time in enumerate(ahc_times[-number_of_steps:]):
        if time <= match_time:
            ahc_times_for_fit_list.append(time)
            coefs_at_different_times_for_fit_list.append(
                coefs_at_different_times[i]
            )
    ahc_times_for_fit = np.array(ahc_times_for_fit_list)
    coefs_at_different_times_for_fit = np.array(
        coefs_at_different_times_for_fit_list
    )
    if ahc_times_for_fit.shape[0] == 0:
        logger.warning(
            "No available AhC times before selected match time; using all"
            " available AhC times, even though numerical errors are likely"
            " larger after the match time"
        )
        ahc_times_for_fit = ahc_times[-number_of_steps:]
        coefs_at_different_times_for_fit = coefs_at_different_times[
            -number_of_steps:
        ]

    logger.info("AhC times available: " + str(ahc_times.shape[0]))
    logger.info(
        "AhC available time range: "
        + str(np.min(ahc_times))
        + " - "
        + str(np.max(ahc_times))
    )
    logger.info("AhC times used: " + str(ahc_times_for_fit.shape[0]))
    logger.info(
        "AhC used time range: "
        + str(np.min(ahc_times_for_fit))
        + " - "
        + str(np.max(ahc_times_for_fit))
    )
    logger.info(
        "Coef times available: " + str(coefs_at_different_times.shape[0])
    )
    logger.info(
        "Coef times used: " + str(coefs_at_different_times_for_fit.shape[0])
    )

    fit_ahc_coefs, fit_ahc_dt_coefs, fit_ahc_dt2_coefs = fit_to_a_cubic(
        ahc_times_for_fit,
        coefs_at_different_times_for_fit,
        match_time,
        zero_coefs,
    )

    # output coefs to H5
    # Note: assumes no translation, so inertial and distorted centers are the
    # same, i.e. both are at the origin. A future update will incorporate
    # translation.

    fit_ahc_coef_mv = ModalVector(fit_ahc_coefs)
    fit_ahc_dt_coef_mv = ModalVector(fit_ahc_dt_coefs)
    fit_ahc_dt2_coef_mv = ModalVector(fit_ahc_dt2_coefs)
    fit_ahc_strahlkorper = Strahlkorper(
        ahc_lmax, ahc_lmax, fit_ahc_coef_mv, ahc_center
    )
    fit_ahc_dt_strahlkorper = Strahlkorper(
        ahc_lmax, ahc_lmax, fit_ahc_dt_coef_mv, ahc_center
    )
    fit_ahc_dt2_strahlkorper = Strahlkorper(
        ahc_lmax, ahc_lmax, fit_ahc_dt2_coef_mv, ahc_center
    )
    legend_ahc, fit_ahc_coefs_to_write = ylm_legend_and_data(
        fit_ahc_strahlkorper, match_time, ahc_lmax
    )
    legend_ahc_dt, fit_ahc_dt_coefs_to_write = ylm_legend_and_data(
        fit_ahc_dt_strahlkorper, match_time, ahc_lmax
    )
    legend_ahc_dt2, fit_ahc_dt2_coefs_to_write = ylm_legend_and_data(
        fit_ahc_dt2_strahlkorper, match_time, ahc_lmax
    )

    for i in range(1, 4, 1):
        legend_ahc[i] = legend_ahc[i].replace("Inertial", "Distorted")
        legend_ahc_dt[i] = legend_ahc_dt[i].replace("Inertial", "Distorted")
        legend_ahc_dt2[i] = legend_ahc_dt2[i].replace("Inertial", "Distorted")

    with spectre_h5.H5File(file_name=path_to_output_h5, mode="a") as h5file:
        ahc_datfile = h5file.insert_dat(
            path="/" + output_subfile_ahc,
            legend=legend_ahc,
            version=0,
        )
        ahc_datfile.append(fit_ahc_coefs_to_write)

    with spectre_h5.H5File(file_name=path_to_output_h5, mode="a") as h5file:
        ahc_dt_datfile = h5file.insert_dat(
            path="/" + output_subfile_dt_ahc,
            legend=legend_ahc_dt,
            version=0,
        )
        ahc_dt_datfile.append(fit_ahc_dt_coefs_to_write)

    with spectre_h5.H5File(file_name=path_to_output_h5, mode="a") as h5file:
        ahc_dt2_datfile = h5file.insert_dat(
            path="/" + output_subfile_dt2_ahc,
            legend=legend_ahc_dt2,
            version=0,
        )
        ahc_dt2_datfile.append(fit_ahc_dt2_coefs_to_write)

    return shape_coefs, info_for_ringdown
