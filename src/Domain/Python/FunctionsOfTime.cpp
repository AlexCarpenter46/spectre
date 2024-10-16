// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Python/FunctionsOfTime.hpp"

#include <memory>
#include <pup.h>
#include <pup_stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <unordered_map>

#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace py = pybind11;

namespace domain::py_bindings {

void bind_functions_of_time(py::module& m) {  // NOLINT
  domain::FunctionsOfTime::register_derived_with_charm();
  py::class_<FunctionsOfTime::FunctionOfTime,
             std::shared_ptr<FunctionsOfTime::FunctionOfTime>>(m,
                                                               "FunctionOfTime")
      .def("time_bounds", &FunctionsOfTime::FunctionOfTime::time_bounds)
      .def("func", &FunctionsOfTime::FunctionOfTime::func, py::arg("t"))
      .def("func_and_deriv", &FunctionsOfTime::FunctionOfTime::func_and_deriv,
           py::arg("t"))
      .def("func_and_2_derivs",
           &FunctionsOfTime::FunctionOfTime::func_and_2_derivs, py::arg("t"));
  m.def(
      "serialize_functions_of_time",
      [](const std::unordered_map<
          std::string,
          std::shared_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
             shared_deserialized_functions_of_time) {
        std::unordered_map<
            std::string,
            std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
            unique_deserialized_functions_of_time{};
        if (shared_deserialized_functions_of_time.contains("Expansion")) {
          unique_deserialized_functions_of_time["Expansion"] =
              shared_deserialized_functions_of_time.at("Expansion")
                  ->get_clone();
        }
        if (shared_deserialized_functions_of_time.contains(
                "ExpansionOuterBoundary")) {
          unique_deserialized_functions_of_time["ExpansionOuterBoundary"] =
              shared_deserialized_functions_of_time
                  .at("ExpansionOuterBoundary")
                  ->get_clone();
        }
        if (shared_deserialized_functions_of_time.contains("Rotation")) {
          unique_deserialized_functions_of_time["Rotation"] =
              shared_deserialized_functions_of_time.at("Rotation")->get_clone();
        }
        if (shared_deserialized_functions_of_time.contains("Translation")) {
          unique_deserialized_functions_of_time["Translation"] =
              shared_deserialized_functions_of_time.at("Translation")
                  ->get_clone();
        }
        return serialize<std::unordered_map<
            std::string,
            std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>>(
            unique_deserialized_functions_of_time);
      },
      py::arg("deserialized_functions_of_time"));
  m.def(
      "deserialize_functions_of_time",
      [](const std::vector<char>& serialized_functions_of_time) {
        return deserialize<std::unordered_map<
            std::string,
            std::shared_ptr<domain::FunctionsOfTime::FunctionOfTime>>>(
            serialized_functions_of_time.data());
      },
      py::arg("serialized_functions_of_time"));
  py::class_<FunctionsOfTime::PiecewisePolynomial<3>,
             FunctionsOfTime::FunctionOfTime,
             std::shared_ptr<FunctionsOfTime::PiecewisePolynomial<3>>>(
      m, "PiecewisePolynomial")
      .def(py::init([](double time,
                       std::array<DataVector, 4>& initial_func_and_derivs,
                       double expiration_time) {
        return FunctionsOfTime::PiecewisePolynomial<3>(
            time, initial_func_and_derivs, expiration_time);
      }))
      .def("time_bounds", &FunctionsOfTime::FunctionOfTime::time_bounds)
      .def("func", &FunctionsOfTime::FunctionOfTime::func, py::arg("t"))
      .def("func_and_deriv", &FunctionsOfTime::FunctionOfTime::func_and_deriv,
           py::arg("t"))
      .def("func_and_2_derivs",
           &FunctionsOfTime::FunctionOfTime::func_and_2_derivs, py::arg("t"));
  py::class_<FunctionsOfTime::QuaternionFunctionOfTime<3>,
             FunctionsOfTime::FunctionOfTime,
             std::shared_ptr<FunctionsOfTime::QuaternionFunctionOfTime<3>>>(
      m, "QuaternionFunctionOfTime")
      .def(
          py::init([](double time, std::array<DataVector, 1>& initial_quat_func,
                      std::array<DataVector, 4>& initial_angle_func,
                      double expiration_time) {
            return FunctionsOfTime::QuaternionFunctionOfTime<3>(
                time, initial_quat_func, initial_angle_func, expiration_time);
          }))
      .def("time_bounds", &FunctionsOfTime::FunctionOfTime::time_bounds)
      .def("func", &FunctionsOfTime::FunctionOfTime::func, py::arg("t"))
      .def("func_and_deriv", &FunctionsOfTime::FunctionOfTime::func_and_deriv,
           py::arg("t"))
      .def("func_and_2_derivs",
           &FunctionsOfTime::FunctionOfTime::func_and_2_derivs, py::arg("t"))
      .def("quat_func", &FunctionsOfTime::FunctionOfTime::func, py::arg("t"))
      .def("quat_func_and_deriv",
           &FunctionsOfTime::FunctionOfTime::func_and_deriv, py::arg("t"))
      .def("quat_func_and_2_derivs",
           &FunctionsOfTime::FunctionOfTime::func_and_2_derivs, py::arg("t"));
}

}  // namespace domain::py_bindings
