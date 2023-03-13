// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/Python/Psi4Real.hpp"

#include <cstddef>
#include <pybind11/pybind11.h>
#include <string>
#include <type_traits>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/GeneralRelativity/Psi4Real.hpp"

namespace py = pybind11;

namespace GeneralRelativity::py_bindings {

namespace {
template <typename Frame>
void bind_psi4real_impl(py::module& m) {  // NOLINT
  m.def(
      "psi4real",
      static_cast<Scalar<DataVector> (*)(const tnsr::ii<DataVector, 3, Frame>&,
                                         const tnsr::ii<DataVector, 3, Frame>&,
                                         const tnsr::ijj<DataVector, 3, Frame>&,
                                         const tnsr::ii<DataVector, 3, Frame>&,
                                         const tnsr::II<DataVector, 3, Frame>&,
                                         const tnsr::I<DataVector, 3, Frame>&)>(
          &::gr::psi_4_real<Frame>),
      py::arg("spatial_ricci"), py::arg("extrinsic_curvature"),
      py::arg("cov_deriv_extrinsic_curvature"), py::arg("spatial_metric"),
      py::arg("inverse_spatial_metric"), py::arg("inertial_coords"));
}
}  // namespace

void bind_psi4real(py::module& m) {  // NOLINT
  bind_psi4real_impl<Frame::Grid>(m);
  bind_psi4real_impl<Frame::Inertial>(m);
}

}  // namespace GeneralRelativity::py_bindings
