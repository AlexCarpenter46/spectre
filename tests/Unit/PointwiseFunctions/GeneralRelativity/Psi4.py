# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np
import math

from .ProjectionOperators import transverse_projection_operator
from .WeylPropagating import weyl_propagating_modes


def psi_4(spatial_ricci, extrinsic_curvature, cov_deriv_extrinsic_curvature,
          spatial_metric, inv_spatial_metric, cartesian_coords):
    raised_r_hat = cartesian_coords / math.sqrt(
        np.einsum("a,b,ab", cartesian_coords, cartesian_coords,
                  spatial_metric))
    r_hat = np.einsum("a,ab", raised_r_hat, spatial_metric)

    inv_projection_tensor = transverse_projection_operator(
        inv_spatial_metric, raised_r_hat)
    projection_tensor = transverse_projection_operator(spatial_metric, r_hat)
    projection_up_lo = np.einsum("ab,ac", inv_projection_tensor,
                                 spatial_metric)

    u8_plus = weyl_propagating_modes(spatial_ricci, extrinsic_curvature,
                                     inv_spatial_metric,
                                     cov_deriv_extrinsic_curvature,
                                     raised_r_hat, inv_projection_tensor,
                                     projection_tensor, projection_up_lo, 1)

    #return np.einsum("a, a", raised_r_hat, r_hat)
    x_coord = np.zeros((3))
    x_coord[0] = cartesian_coords[0]
    x_hat = x_coord / math.sqrt(
        np.einsum("a,b,ab", x_coord, x_coord, spatial_metric))
    y_coord = np.zeros((3))
    y_coord[1] = cartesian_coords[1]
    y_hat = y_coord / math.sqrt(
        np.einsum("a,b,ab", y_coord, y_coord, spatial_metric))

    return (0.5 * np.einsum("ab, ab", u8_plus,
                            ((np.einsum("a,b", y_hat, y_hat)) -
                             (np.einsum("a,b", x_hat, x_hat)))))
