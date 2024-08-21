// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/Python/Rectilinear.hpp"

#include <array>
#include <cstddef>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/Rectilinear.hpp"

namespace py = pybind11;

namespace domain::creators::py_bindings {
void bind_brick(py::module& m) {
  py::class_<Brick, DomainCreator<3>>(m, "Brick")
      .def(py::init<std::array<double, 3>, std::array<double, 3>,
                    std::array<size_t, 3>, std::array<size_t, 3>,
                    std::array<bool, 3>>(),
           py::arg("lower_xyz"), py::arg("upper_xyz"),
           py::arg("initial_refinement_level_xyz"),
           py::arg("initial_number_of_grid_points_in_xyz"),
           py::arg("is_periodic_in_xyz"));
}
}  // namespace domain::creators::py_bindings
