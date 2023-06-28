// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"

#include <cmath>
#include <ostream>
#include <pup.h>
#include <pup_stl.h>
#include <type_traits>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Expressions/Evaluate.hpp"
#include "DataStructures/Tensor/Identity.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "NumericalAlgorithms/RootFinding/RootBracketing.hpp"
#include "NumericalAlgorithms/RootFinding/TOMS748.hpp"
#include "PointwiseFunctions/MathFunctions/Gaussian.hpp"
#include "PointwiseFunctions/MathFunctions/MathFunction.hpp"
#include "Utilities/DereferenceWrapper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/StdHelpers.hpp"

namespace domain::CoordinateMaps::TimeDependent {

template <size_t Dim>
Translation<Dim>::Translation(
    std::string function_of_time_name,
    std::unique_ptr<MathFunction<1, Frame::Inertial>>
        radial_function,
    std::array<double, Dim>& center)
    : f_of_t_name_(std::move(function_of_time_name)),
      f_of_r_(std::move(radial_function)),
      center_(center) {}

template <size_t Dim>
template <typename DataType>
std::array<tt::remove_cvref_wrap_t<DataType>, Dim> Translation<Dim>::operator()(
    const std::array<DataType, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.count(f_of_t_name_) == 1,
         "The function of time '" << f_of_t_name_
                                  << "' is not one of the known functions of "
                                     "time. The known functions of time are: "
                                  << keys_of(functions_of_time));
  std::array<tt::remove_cvref_wrap_t<DataType>, Dim> result{};
  const DataVector function_of_time =
      functions_of_time.at(f_of_t_name_)->func(time)[0];
  ASSERT(function_of_time.size() == Dim,
         "The dimension of the function of time ("
             << function_of_time.size()
             << ") does not match the dimension of the translation map (" << Dim
             << ").");
  // finding the r coord
  std::array<tt::remove_cvref_wrap_t<DataType>, Dim> distance_to_center{};
  // for (size_t i = 0; i < Dim; i++) {
  //     distance_to_center[i] = source_coords[i];
  // }
  Scalar<tt::remove_cvref_wrap_t<DataType>> radius(
      get_size(distance_to_center[0]));
  for (size_t i = 0; i < Dim; i++) {
    distance_to_center[i] = source_coords[i] - center_[i];
    radius.get() += square(distance_to_center[i]);
  }
  radius.get() = sqrt(radius.get());
  for (size_t i = 0; i < Dim; i++) {
    gsl::at(result, i) = gsl::at(source_coords, i) +
                         function_of_time[i] * ((*f_of_r_)(radius.get()));
  }
  return result;
}

template <size_t Dim>
std::optional<std::array<double, Dim>> Translation<Dim>::inverse(
    const std::array<double, Dim>& translated_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.count(f_of_t_name_) == 1,
         "The function of time '" << f_of_t_name_
                                  << "' is not one of the known functions of "
                                     "time. The known functions of time are: "
                                  << keys_of(functions_of_time));
  std::array<double, Dim> result{};
  const DataVector function_of_time =
      functions_of_time.at(f_of_t_name_)->func(time)[0];
  ASSERT(function_of_time.size() == Dim,
         "The dimension of the function of time ("
             << function_of_time.size()
             << ") does not match the dimension of the translation map (" << Dim
             << ").");
  // finding the shifted r coord and square of the translation map
  // that'll be used for the root finder.
  std::array<double, Dim> distance_to_center = translated_coords;
  // change this to just a double cause we don't need a Scalar here.
  Scalar<double> magnitude_function_of_time{0.};
  double radius_squared = 0.;
  for (size_t i = 0; i < Dim; i++) {
    // shifted center
    distance_to_center[i] -= center_[i];
    radius_squared += square(distance_to_center[i]);
    // square of the function of time
    magnitude_function_of_time.get() += square(function_of_time[i]);
  }
  const double translated_radius = sqrt(radius_squared);
  const double& root = root_finder(translated_coords, translated_radius,
                                   magnitude_function_of_time);
  for (size_t i = 0; i < Dim; i++) {
    gsl::at(result, i) =
        gsl::at(translated_coords, i) - function_of_time[i] * (*f_of_r_)(root);
  }
  return result;
}

template <size_t Dim>
template <typename DataType>
std::array<tt::remove_cvref_wrap_t<DataType>, Dim>
Translation<Dim>::frame_velocity(
    const std::array<DataType, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.count(f_of_t_name_) == 1,
         "The function of time '" << f_of_t_name_
                                  << "' is not one of the known functions of "
                                     "time. The known functions of time are: "
                                  << keys_of(functions_of_time));
  std::array<tt::remove_cvref_wrap_t<DataType>, Dim> result{};
  const auto function_of_time_and_deriv =
      functions_of_time.at(f_of_t_name_)->func_and_deriv(time);
  const DataVector& velocity = function_of_time_and_deriv[1];
  ASSERT(velocity.size() == Dim,
         "The dimension of the function of time ("
             << velocity.size()
             << ") does not match the dimension of the translation map (" << Dim
             << ").");
  for (size_t i = 0; i < Dim; i++) {
    gsl::at(result, i) = make_with_value<tt::remove_cvref_wrap_t<DataType>>(
        dereference_wrapper(gsl::at(source_coords, i)), velocity[i]);
  }
  return result;
}

template <size_t Dim>
template <typename DataType>
tnsr::Ij<tt::remove_cvref_wrap_t<DataType>, Dim, Frame::NoFrame>
Translation<Dim>::jacobian(
    const std::array<DataType, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  std::array<tt::remove_cvref_wrap_t<DataType>, Dim> distance_to_center{};
  Scalar<tt::remove_cvref_wrap_t<DataType>> radius(
      get_size(distance_to_center[0]));
  for (size_t i = 0; i < Dim; i++) {
    distance_to_center[i] = source_coords[i] - center_[i];
    radius.get() += square(distance_to_center[i]);
  }
  radius.get() = sqrt(radius.get());
  // radius = sqrt(radius);

  const DataVector function_of_time =
      functions_of_time.at(f_of_t_name_)->func(time)[0];

  auto result = make_with_value<
      tnsr::Ij<tt::remove_cvref_wrap_t<DataType>, Dim, Frame::NoFrame>>(
      dereference_wrapper(source_coords[0]), 0.0);
  for (size_t i = 0; i < Dim; i++) {
    for (size_t j = 0; j < Dim; j++) {
      if constexpr (std::is_same_v<tt::remove_cvref_wrap_t<DataType>, double>) {
        if (radius.get() > 1.e-13) {
          result.get(i, j) = (*f_of_r_).first_deriv(radius.get()) *
                             gsl::at(function_of_time, i) *
                             gsl::at(distance_to_center, j) / radius.get();
        } else {
          result.get(i, j) = (*f_of_r_).second_deriv(radius.get()) *
                             gsl::at(function_of_time, i) *
                             gsl::at(distance_to_center, j);
        }
      } else {
        for (size_t k = 0; k < radius.size(); k++) {
          if (radius.get()[k] > 1.e-13) {
            result.get(i, j)[k] = (*f_of_r_).first_deriv(radius.get()[k]) *
                                  gsl::at(function_of_time, i) *
                                  gsl::at(distance_to_center, j)[k] /
                                  radius.get()[k];
          } else {
            result.get(i, j)[k] = (*f_of_r_).second_deriv(radius.get()[k]) *
                                  gsl::at(function_of_time, i) *
                                  gsl::at(distance_to_center, j)[k];
          }
        }
      }
    }
    result.get(i, i) += 1.0;
  }
  return result;
}

template <size_t Dim>
template <typename DataType>
tnsr::Ij<tt::remove_cvref_wrap_t<DataType>, Dim, Frame::NoFrame>
Translation<Dim>::inv_jacobian(
    const std::array<DataType, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  return determinant_and_inverse(
             jacobian(source_coords, time, functions_of_time))
      .second;
}

template <size_t Dim>
double Translation<Dim>::root_finder(
    const std::array<double, Dim>& translated_coords,
    const double& translated_radius,
    const Scalar<double>& magnitude_function_of_time) const {
  const double center_offset = sqrt(magnitude_function_of_time.get());
  double lower_bound = (translated_radius - center_offset) * (1.0 - 1.e-9);
  double upper_bound = (translated_radius + center_offset) * (1.0 + 1.e-9);
  double function_at_lower_bound = (*f_of_r_)(lower_bound);
  double function_at_upper_bound = (*f_of_r_)(upper_bound);
  // RootFinder::bracket_possibly_undefined_function_in_interval(
  //     &lower_bound, &upper_bound, &function_at_lower_bound,
  //     &function_at_upper_bound, *f_of_r_);
  const double absolute_tol = std::numeric_limits<double>::epsilon() *
                              std::max(translated_radius, center_offset);
  const double relative_tol = std::numeric_limits<double>::epsilon();
  return RootFinder::toms748(*f_of_r_, lower_bound, upper_bound,
                             function_at_lower_bound, function_at_upper_bound,
                             absolute_tol, relative_tol);
}

template <size_t Dim>
void Translation<Dim>::pup(PUP::er& p) {
  size_t version = 1;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | f_of_t_name_;
  }
  if (version >= 1) {
    p | f_of_r_;
    p | center_;
  }
}

template <size_t Dim>
bool operator==(const Translation<Dim>& lhs, const Translation<Dim>& rhs) {
  return lhs.f_of_t_name_ == rhs.f_of_t_name_;
}

// Explicit instantiations
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                              \
  template class Translation<DIM(data)>;                  \
  template bool operator==(const Translation<DIM(data)>&, \
                           const Translation<DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data)                                           \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)> \
  Translation<DIM(data)>::operator()(                                  \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;                                    \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)> \
  Translation<DIM(data)>::frame_velocity(                              \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;                                    \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),   \
                    Frame::NoFrame>                                    \
  Translation<DIM(data)>::jacobian(                                    \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;                                    \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),   \
                    Frame::NoFrame>                                    \
  Translation<DIM(data)>::inv_jacobian(                                \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3),
                        (double, DataVector,
                         std::reference_wrapper<const double>,
                         std::reference_wrapper<const DataVector>))
#undef DIM
#undef DTYPE
#undef INSTANTIATE
}  // namespace domain::CoordinateMaps::TimeDependent
