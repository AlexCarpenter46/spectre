# Distributed under the MIT License.
# See LICENSE.txt for details.

from spectre.GeneralRelativity import psi4real
from spectre.DataStructures import DataVector
from spectre.DataStructures.Tensor import tnsr, Frame
from spectre.Spectral import Basis, Mesh, Quadrature, collocation_points

import numpy as np
import unittest


class TestPsi4(unittest.TestCase):
    def test_psi4(self):
        mesh = Mesh[1](4, Basis.Legendre, Quadrature.GaussLobatto)
        x = collocation_points(mesh)
        inertial_coordinates = tnsr.I[DataVector, 1, Frame.Inertial]([x])
        spatial_ricci = tnsr.ij[DataVector, 3, Frame.Inertial](num_point=9,
                                                               fill=2.)


if __name__ == '__main__':
    unittest.main()
