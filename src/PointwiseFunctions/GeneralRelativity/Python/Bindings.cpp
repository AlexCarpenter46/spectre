// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/pybind11.h>

#include "PointwiseFunctions/GeneralRelativity/Python/Psi4Real.hpp"

namespace py = pybind11;

namespace GeneralRelativity {

PYBIND11_MODULE(_PyGeneralRelativity, m) {  // NOLINT
  py::module_::import("spectre.DataStructures");
  py::module_::import("spectre.DataStructures.Tensor");
  py::module_::import("spectre.Spectral");
  py_bindings::bind_psi4real(m);
}
}  // namespace GeneralRelativity
