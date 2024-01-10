// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Identity.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.tpp"
#include "Domain/CoordinateMaps/TimeDependent/RotScaleTrans.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/CoordinateMaps/TimeDependent/RotationMatrixHelpers.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/Domain/CoordinateMaps/TestMapHelpers.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GetOutput.hpp"

template <size_t Dim>
void test_RotScaleTrans() {
  MAKE_GENERATOR(gen);
  // define vars for FunctionOfTime::PiecewisePolynomial f(t) = t**2.
  double t = -1.0;
  const double dt = 0.6;
  const double final_time = 4.0;
  constexpr size_t deriv_order = 3;
  const double inner_radius = 1.0;
  const double outer_radius = 50.0;
  const double angle = 1.0;
  const double omega = -2.0;
  const double dtomega = 2.0;
  const double d2tomega = 0.0;

  using Polynomial = domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>;
  using FoftPtr = std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>;
  using QuatFoT =
      domain::FunctionsOfTime::QuaternionFunctionOfTime<deriv_order>;
  std::unordered_map<std::string, FoftPtr> f_of_t_list{};
  std::array<DataVector, deriv_order + 1> init_func{};
  const std::array<DataVector, deriv_order + 1> init_func_trans{
      {{Dim, 1.0}, {Dim, -2.0}, {Dim, 2.0}, {Dim, 0.0}}};
  const std::array<DataVector, deriv_order + 1> init_func_a{
      {{0.96}, {0.0}, {0.0}, {0.0}}};
  const std::array<DataVector, deriv_order + 1> init_func_b{
      {{0.58}, {0.0}, {0.0}, {0.0}}};
  const std::string rot_f_of_t_name{"rotation_angle"};
  Approx custom_approx = Approx::custom().epsilon(1.0).scale(1.0);

  if constexpr (Dim == 2) {
    init_func = {{{angle}, {omega}, {dtomega}, {d2tomega}}};
    f_of_t_list[rot_f_of_t_name] =
        std::make_unique<Polynomial>(t, init_func, final_time + dt);
    custom_approx = Approx::custom().epsilon(5.0e-13).scale(1.0);
  } else {
    // Axis of rotation nhat = (1.0, -1.0, 1.0) / sqrt(3.0)
    DataVector axis{{1.0, -1.0, 1.0}};
    axis /= sqrt(3.0);
    init_func = {axis * angle, axis * omega, axis * dtomega, axis * d2tomega};
    // initial quaternion is (cos(angle/2), nhat*sin(angle/2))
    const std::array<DataVector, 1> init_quat{
        DataVector{{cos(angle / 2.0), axis[0] * sin(angle / 2.0),
                    axis[1] * sin(angle / 2.0), axis[2] * sin(angle / 2.0)}}};
    f_of_t_list[rot_f_of_t_name] =
        std::make_unique<QuatFoT>(t, init_quat, init_func, final_time + dt);
    custom_approx = Approx::custom().epsilon(5.0e-11).scale(1.0);
  }
  f_of_t_list["expansion_a"] =
      std::make_unique<Polynomial>(t, init_func_a, final_time);
  f_of_t_list["expansion_b"] =
      std::make_unique<Polynomial>(t, init_func_b, final_time);
  f_of_t_list["translation"] =
      std::make_unique<Polynomial>(t, init_func_trans, final_time + dt);

  const FoftPtr& scale_a_of_t = f_of_t_list.at("expansion_a");
  const FoftPtr& scale_b_of_t = f_of_t_list.at("expansion_b");
  const FoftPtr& trans_f_of_t = f_of_t_list.at("translation");
  const std::pair<std::string, std::string> scale_pair{"expansion_a",
                                                       "expansion_b"};
  // Rotation, Scaling, Translation Non-Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      rot_scale_trans_map_non_rigid{scale_pair,    "rotation_angle",
                                    "translation", inner_radius,
                                    outer_radius,  false};

  // Rotation, Scaling, Translation Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      rot_scale_trans_map_rigid{scale_pair,   "rotation_angle", "translation",
                                inner_radius, outer_radius,     true};

  // Rotation, Scaling Non-Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      rot_scale_map_non_rigid{scale_pair,   "rotation_angle", std::nullopt,
                              inner_radius, outer_radius,     false};

  // Rotation, Translation Non-Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      rot_trans_map_non_rigid{std::nullopt, "rotation_angle", "translation",
                              inner_radius, outer_radius,     false};

  // Rotation, Translation Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      rot_trans_map_rigid{std::nullopt, "rotation_angle", "translation",
                          inner_radius, outer_radius,     true};

  // Scaling, Translation Non-Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      scale_trans_map_non_rigid{scale_pair,   std::nullopt, "translation",
                                inner_radius, outer_radius, false};

  // Scaling, Translation Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      scale_trans_map_rigid{scale_pair,   std::nullopt, "translation",
                            inner_radius, outer_radius, true};

  // Rotation
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim> rot_map{
      std::nullopt, "rotation_angle", std::nullopt,
      inner_radius, outer_radius,     true};

  // Scaling
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      scale_map_non_rigid{scale_pair,   std::nullopt, std::nullopt,
                          inner_radius, outer_radius, false};

  // Translation Non-Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      trans_map_non_rigid{std::nullopt, std::nullopt, "translation",
                          inner_radius, outer_radius, false};

  // Translation Rigid
  const domain::CoordinateMaps::TimeDependent::RotScaleTrans<Dim>
      trans_map_rigid{std::nullopt, std::nullopt, "translation",
                      inner_radius, outer_radius, true};

  UniformCustomDistribution<double> dist_double{-5.0, 5.0};
  UniformCustomDistribution<double> far_dist_double{50.0, 60.0};
  std::array<double, Dim> point_xi{};
  std::array<DataVector, Dim> point_xi_dv{};
  std::array<double, Dim> far_point_xi{};
  fill_with_random_values(make_not_null(&point_xi), make_not_null(&gen),
                          make_not_null(&dist_double));
  fill_with_random_values(make_not_null(&far_point_xi), make_not_null(&gen),
                          make_not_null(&far_dist_double));
  for (size_t i = 0; i < Dim; i++) {
    auto points = make_with_random_values<DataVector>(
        make_not_null(&gen), make_not_null(&dist_double), DataVector(5));
    gsl::at(point_xi_dv, i) = points;
  }

  while (t < final_time) {
    std::array<double, Dim> translation{};
    std::array<double, Dim> deriv_translation{};
    std::array<double, Dim> expected_rotation{};
    std::array<double, Dim> expected_rotation_deriv{};
    std::array<double, Dim> far_expected_rotation{};
    std::array<double, Dim> far_expected_rotation_deriv{};
    const double radius = magnitude(point_xi);
    const DataVector radius_dv = magnitude(point_xi_dv);
    std::array<double, Dim> frame_vel_rot_scale{};
    std::array<double, Dim> far_frame_vel_rot_scale{};
    std::array<double, Dim> frame_vel_rot_scale_rigid{};
    std::array<double, Dim> frame_vel_t_nr{0.};
    const Matrix rot_matrix =
        rotation_matrix<Dim>(t, *(f_of_t_list[rot_f_of_t_name]));
    const Matrix deriv_rot_matrix =
        rotation_matrix_deriv<Dim>(t, *(f_of_t_list[rot_f_of_t_name]));
    const double scale_a = scale_a_of_t->func_and_deriv(t)[0][0];
    const double scale_b = scale_b_of_t->func_and_deriv(t)[0][0];
    const double d_scale_a = scale_a_of_t->func_and_deriv(t)[1][0];
    const double d_scale_b = scale_b_of_t->func_and_deriv(t)[1][0];
    double radial_scaling_factor = 0.0;
    double deriv_radial_scaling_factor = 0.0;
    double radial_translation_factor = 0.0;
    custom_approx = Approx::custom().epsilon(1.e-10);
    for (size_t i = 0; i < Dim; i++) {
      gsl::at(translation, i) = square(t);
      gsl::at(deriv_translation, i) = trans_f_of_t->func_and_deriv(t)[1][i];
      for (size_t j = 0; j < Dim; j++) {
        gsl::at(expected_rotation, i) +=
            rot_matrix(i, j) * gsl::at(point_xi, j);
        gsl::at(expected_rotation_deriv, i) +=
            deriv_rot_matrix(i, j) * gsl::at(point_xi, j);
        gsl::at(frame_vel_rot_scale_rigid, i) +=
            gsl::at(point_xi, j) *
            (scale_a * deriv_rot_matrix(i, j) + d_scale_a * rot_matrix(i, j));
        gsl::at(far_expected_rotation, i) +=
            rot_matrix(i, j) * gsl::at(far_point_xi, j);
        gsl::at(far_expected_rotation_deriv, i) +=
            deriv_rot_matrix(i, j) * gsl::at(far_point_xi, j);
        gsl::at(far_frame_vel_rot_scale, i) +=
            gsl::at(far_point_xi, j) *
            (scale_b * deriv_rot_matrix(i, j) + d_scale_b * rot_matrix(i, j));
      }
    }
    // Operator
    CHECK_ITERABLE_APPROX(trans_map_rigid(point_xi, t, f_of_t_list),
                          point_xi + translation);
    CHECK_ITERABLE_APPROX(rot_map(point_xi, t, f_of_t_list), expected_rotation);
    CHECK_ITERABLE_APPROX(rot_trans_map_rigid(point_xi, t, f_of_t_list),
                          expected_rotation + translation);
    CHECK_ITERABLE_APPROX(scale_trans_map_rigid(point_xi, t, f_of_t_list),
                          point_xi * scale_a + translation);
    CHECK_ITERABLE_APPROX(rot_scale_trans_map_rigid(point_xi, t, f_of_t_list),
                          expected_rotation * scale_a + translation);
    // Frame Velocity
    CHECK_ITERABLE_APPROX(rot_map.frame_velocity(point_xi, t, f_of_t_list),
                          expected_rotation_deriv);
    CHECK_ITERABLE_APPROX(
        trans_map_rigid.frame_velocity(point_xi, t, f_of_t_list),
        deriv_translation);
    CHECK_ITERABLE_APPROX(
        rot_trans_map_rigid.frame_velocity(point_xi, t, f_of_t_list),
        expected_rotation_deriv + deriv_translation);
    CHECK_ITERABLE_APPROX(
        scale_trans_map_rigid.frame_velocity(point_xi, t, f_of_t_list),
        point_xi * (scale_a + d_scale_a) + deriv_translation);
    // Inverse
    CHECK_ITERABLE_APPROX(
        rot_map.inverse(expected_rotation, t, f_of_t_list).value(), point_xi);
    CHECK_ITERABLE_APPROX(
        trans_map_rigid.inverse(point_xi + translation, t, f_of_t_list).value(),
        point_xi);
    CHECK_ITERABLE_APPROX(
        rot_trans_map_rigid
            .inverse(expected_rotation + translation, t, f_of_t_list)
            .value(),
        point_xi);
    CHECK_ITERABLE_APPROX(
        scale_trans_map_rigid
            .inverse(point_xi * scale_a + translation, t, f_of_t_list)
            .value(),
        point_xi);
    CHECK_ITERABLE_APPROX(
        rot_scale_trans_map_rigid
            .inverse(expected_rotation * scale_a + translation, t, f_of_t_list)
            .value(),
        point_xi);

    if (radius <= inner_radius) {
      frame_vel_t_nr = deriv_translation;
      frame_vel_rot_scale = frame_vel_rot_scale_rigid;
      CHECK_ITERABLE_APPROX(scale_map_non_rigid(point_xi, t, f_of_t_list),
                            point_xi * scale_a);
      CHECK_ITERABLE_APPROX(trans_map_non_rigid(point_xi, t, f_of_t_list),
                            point_xi + translation);
      CHECK_ITERABLE_APPROX(rot_scale_map_non_rigid(point_xi, t, f_of_t_list),
                            expected_rotation * scale_a);
      CHECK_ITERABLE_APPROX(rot_trans_map_non_rigid(point_xi, t, f_of_t_list),
                            expected_rotation + translation);
      CHECK_ITERABLE_APPROX(scale_trans_map_non_rigid(point_xi, t, f_of_t_list),
                            point_xi * scale_a + translation);
      CHECK_ITERABLE_APPROX(
          rot_scale_trans_map_non_rigid(point_xi, t, f_of_t_list),
          expected_rotation * scale_a + translation);
      // Frame Velocity
      CHECK_ITERABLE_APPROX(
          scale_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
          point_xi * (scale_a + d_scale_a));
      CHECK_ITERABLE_APPROX(
          trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
          deriv_translation);
      CHECK_ITERABLE_APPROX(
          rot_trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
          expected_rotation_deriv + deriv_translation);
      CHECK_ITERABLE_APPROX(
          scale_trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
          point_xi * (scale_a + d_scale_a) + deriv_translation);
      // Inverse
      CHECK_ITERABLE_APPROX(
          scale_map_non_rigid.inverse(point_xi * scale_a, t, f_of_t_list)
              .value(),
          point_xi);
      CHECK_ITERABLE_APPROX(
          trans_map_non_rigid.inverse(point_xi + translation, t, f_of_t_list)
              .value(),
          point_xi);
      CHECK_ITERABLE_APPROX(
          rot_scale_map_non_rigid
              .inverse(expected_rotation * scale_a, t, f_of_t_list)
              .value(),
          point_xi);
      CHECK_ITERABLE_APPROX(
          rot_trans_map_non_rigid
              .inverse(expected_rotation + translation, t, f_of_t_list)
              .value(),
          point_xi);
      CHECK_ITERABLE_APPROX(
          scale_trans_map_non_rigid
              .inverse(point_xi * scale_a + translation, t, f_of_t_list)
              .value(),
          point_xi);
      CHECK_ITERABLE_APPROX(
          rot_scale_trans_map_non_rigid
              .inverse(expected_rotation * scale_a + translation, t,
                       f_of_t_list)
              .value(),
          point_xi);
    } else if (radius > inner_radius and radius < outer_radius) {
      if (1 - radius / outer_radius < .1) {
        radial_scaling_factor =
            ((outer_radius - radius) * (scale_a - scale_b) * inner_radius) /
            ((outer_radius - inner_radius) * radius);
        deriv_radial_scaling_factor =
            ((outer_radius - radius) * (d_scale_a - d_scale_b) * inner_radius) /
            ((outer_radius - inner_radius) * radius);
        radial_translation_factor =
            (outer_radius - radius) / (outer_radius - inner_radius);
        frame_vel_t_nr = deriv_translation * radial_translation_factor;
        for (size_t i = 0; i < Dim; i++) {
          for (size_t j = 0; j < Dim; j++) {
            gsl::at(frame_vel_rot_scale, i) +=
                gsl::at(point_xi, j) *
                (deriv_rot_matrix(i, j) * (scale_b + radial_scaling_factor) +
                 rot_matrix(i, j) * (d_scale_b + deriv_radial_scaling_factor));
          }
        }
        CHECK_ITERABLE_APPROX(scale_map_non_rigid(point_xi, t, f_of_t_list),
                              point_xi * (radial_scaling_factor + scale_b));
        CHECK_ITERABLE_APPROX(
            trans_map_non_rigid(point_xi, t, f_of_t_list),
            point_xi + translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            rot_scale_map_non_rigid(point_xi, t, f_of_t_list),
            expected_rotation * (radial_scaling_factor + scale_b));
        CHECK_ITERABLE_APPROX(
            rot_trans_map_non_rigid(point_xi, t, f_of_t_list),
            expected_rotation + translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            scale_trans_map_non_rigid(point_xi, t, f_of_t_list),
            point_xi * (radial_scaling_factor + scale_b) +
                translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            rot_scale_trans_map_non_rigid(point_xi, t, f_of_t_list),
            expected_rotation * (radial_scaling_factor + scale_b) +
                translation * radial_translation_factor);
        // Frame Velocity
        CHECK_ITERABLE_APPROX(
            scale_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            point_xi * (scale_b + d_scale_b + radial_scaling_factor +
                        deriv_radial_scaling_factor));
        CHECK_ITERABLE_APPROX(
            trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            deriv_translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            rot_trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            expected_rotation_deriv +
                deriv_translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            scale_trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            point_xi * (scale_b + d_scale_b + radial_scaling_factor +
                        deriv_radial_scaling_factor) +
                deriv_translation * radial_translation_factor);
        // Inverse
        CHECK_ITERABLE_CUSTOM_APPROX(
            scale_map_non_rigid
                .inverse(point_xi * (radial_scaling_factor + scale_b), t,
                         f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            trans_map_non_rigid
                .inverse(point_xi + translation * radial_translation_factor, t,
                         f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            rot_scale_map_non_rigid
                .inverse(expected_rotation * (radial_scaling_factor + scale_b),
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            rot_trans_map_non_rigid
                .inverse(
                    expected_rotation + translation * radial_translation_factor,
                    t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            scale_trans_map_non_rigid
                .inverse(point_xi * (radial_scaling_factor + scale_b) +
                             translation * radial_translation_factor,
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            rot_scale_trans_map_non_rigid
                .inverse(expected_rotation * (radial_scaling_factor + scale_b) +
                             translation * radial_translation_factor,
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
      } else {
        radial_scaling_factor =
            ((inner_radius - radius) * (scale_a - scale_b) * outer_radius) /
            ((outer_radius - inner_radius) * radius);
        deriv_radial_scaling_factor =
            ((inner_radius - radius) * (d_scale_a - d_scale_b) * outer_radius) /
            ((outer_radius - inner_radius) * radius);
        radial_translation_factor =
            (inner_radius - radius) / (outer_radius - inner_radius);
        frame_vel_t_nr =
            deriv_translation + deriv_translation * radial_translation_factor;
        for (size_t i = 0; i < Dim; i++) {
          for (size_t j = 0; j < Dim; j++) {
            gsl::at(frame_vel_rot_scale, i) +=
                gsl::at(point_xi, j) *
                (deriv_rot_matrix(i, j) * (scale_a + radial_scaling_factor) +
                 rot_matrix(i, j) * (d_scale_a + deriv_radial_scaling_factor));
          }
        }
        CHECK_ITERABLE_APPROX(scale_map_non_rigid(point_xi, t, f_of_t_list),
                              point_xi * (radial_scaling_factor + scale_a));
        CHECK_ITERABLE_APPROX(
            trans_map_non_rigid(point_xi, t, f_of_t_list),
            point_xi + translation + translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            rot_scale_map_non_rigid(point_xi, t, f_of_t_list),
            expected_rotation * (radial_scaling_factor + scale_a));
        CHECK_ITERABLE_APPROX(rot_trans_map_non_rigid(point_xi, t, f_of_t_list),
                              expected_rotation + translation +
                                  translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            scale_trans_map_non_rigid(point_xi, t, f_of_t_list),
            point_xi * (radial_scaling_factor + scale_a) + translation +
                translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            rot_scale_trans_map_non_rigid(point_xi, t, f_of_t_list),
            expected_rotation * (radial_scaling_factor + scale_a) +
                translation + translation * radial_translation_factor);
        // Frame Velocity
        CHECK_ITERABLE_APPROX(
            scale_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            point_xi * (scale_a + d_scale_a + radial_scaling_factor +
                        deriv_radial_scaling_factor));
        CHECK_ITERABLE_APPROX(
            trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            deriv_translation + deriv_translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            rot_trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            expected_rotation_deriv + deriv_translation +
                deriv_translation * radial_translation_factor);
        CHECK_ITERABLE_APPROX(
            scale_trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
            point_xi * (scale_a + d_scale_a + radial_scaling_factor +
                        deriv_radial_scaling_factor) +
                deriv_translation +
                deriv_translation * radial_translation_factor);
        // Inverse
        CHECK_ITERABLE_CUSTOM_APPROX(
            scale_map_non_rigid
                .inverse(point_xi * (radial_scaling_factor + scale_a), t,
                         f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            trans_map_non_rigid
                .inverse(point_xi + translation +
                             translation * radial_translation_factor,
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            rot_scale_map_non_rigid
                .inverse(expected_rotation * (radial_scaling_factor + scale_a),
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            rot_trans_map_non_rigid
                .inverse(expected_rotation + translation +
                             translation * radial_translation_factor,
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            scale_trans_map_non_rigid
                .inverse(point_xi * (radial_scaling_factor + scale_a) +
                             translation +
                             translation * radial_translation_factor,
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            rot_scale_trans_map_non_rigid
                .inverse(expected_rotation * (radial_scaling_factor + scale_a) +
                             translation +
                             translation * radial_translation_factor,
                         t, f_of_t_list)
                .value(),
            point_xi, custom_approx);
      }
    }
    // Section for testing points on/beyond outer boundary.
    CHECK_ITERABLE_APPROX(scale_map_non_rigid(far_point_xi, t, f_of_t_list),
                          far_point_xi * scale_b);
    CHECK_ITERABLE_APPROX(trans_map_non_rigid(far_point_xi, t, f_of_t_list),
                          far_point_xi);
    CHECK_ITERABLE_APPROX(rot_scale_map_non_rigid(far_point_xi, t, f_of_t_list),
                          far_expected_rotation * scale_b);
    CHECK_ITERABLE_APPROX(rot_trans_map_non_rigid(far_point_xi, t, f_of_t_list),
                          far_expected_rotation);
    CHECK_ITERABLE_APPROX(
        rot_scale_trans_map_non_rigid(far_point_xi, t, f_of_t_list),
        far_expected_rotation * scale_b);
    // Frame Velocity
    CHECK_ITERABLE_APPROX(
        scale_map_non_rigid.frame_velocity(far_point_xi, t, f_of_t_list),
        far_point_xi * (scale_b + d_scale_b));
    CHECK_ITERABLE_APPROX(
        rot_trans_map_non_rigid.frame_velocity(far_point_xi, t, f_of_t_list),
        far_expected_rotation_deriv);
    CHECK_ITERABLE_APPROX(
        scale_trans_map_non_rigid.frame_velocity(far_point_xi, t, f_of_t_list),
        far_point_xi * (scale_b + d_scale_b));
    // Inverse
    CHECK_ITERABLE_APPROX(
        scale_map_non_rigid.inverse(far_point_xi * scale_b, t, f_of_t_list)
            .value(),
        far_point_xi);
    CHECK_ITERABLE_APPROX(
        trans_map_non_rigid.inverse(far_point_xi, t, f_of_t_list).value(),
        far_point_xi);
    CHECK_ITERABLE_APPROX(
        rot_scale_map_non_rigid
            .inverse(far_expected_rotation * scale_b, t, f_of_t_list)
            .value(),
        far_point_xi);
    CHECK_ITERABLE_APPROX(
        rot_trans_map_non_rigid.inverse(far_expected_rotation, t, f_of_t_list)
            .value(),
        far_point_xi);
    CHECK_ITERABLE_APPROX(scale_trans_map_non_rigid
                              .inverse(far_point_xi * scale_b, t, f_of_t_list)
                              .value(),
                          far_point_xi);
    CHECK_ITERABLE_APPROX(
        rot_scale_trans_map_non_rigid
            .inverse(far_expected_rotation * scale_b, t, f_of_t_list)
            .value(),
        far_point_xi);

    // frame velocity from different sections
    CHECK_ITERABLE_APPROX(
        rot_scale_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
        frame_vel_rot_scale);
    CHECK_ITERABLE_APPROX(
        rot_scale_trans_map_rigid.frame_velocity(point_xi, t, f_of_t_list),
        frame_vel_rot_scale_rigid + deriv_translation);
    CHECK_ITERABLE_APPROX(
        rot_scale_trans_map_non_rigid.frame_velocity(point_xi, t, f_of_t_list),
        frame_vel_rot_scale + frame_vel_t_nr);

    if (radius <= inner_radius * .99 or radius >= inner_radius * 1.01) {
      test_jacobian(rot_map, point_xi, t, f_of_t_list);
      test_inv_jacobian(rot_map, point_xi, t, f_of_t_list);
      test_jacobian(scale_map_non_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(scale_map_non_rigid, point_xi, t, f_of_t_list);
      test_jacobian(trans_map_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(trans_map_rigid, point_xi, t, f_of_t_list);
      test_jacobian(trans_map_non_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(trans_map_non_rigid, point_xi, t, f_of_t_list);
      test_jacobian(rot_scale_map_non_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(rot_scale_map_non_rigid, point_xi, t, f_of_t_list);
      test_jacobian(rot_trans_map_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(rot_trans_map_rigid, point_xi, t, f_of_t_list);
      test_jacobian(rot_trans_map_non_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(rot_trans_map_non_rigid, point_xi, t, f_of_t_list);
      test_jacobian(scale_trans_map_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(scale_trans_map_rigid, point_xi, t, f_of_t_list);
      test_jacobian(scale_trans_map_non_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(scale_trans_map_non_rigid, point_xi, t, f_of_t_list);
      test_jacobian(rot_scale_trans_map_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(rot_scale_trans_map_rigid, point_xi, t, f_of_t_list);
      test_jacobian(rot_scale_trans_map_non_rigid, point_xi, t, f_of_t_list);
      test_inv_jacobian(rot_scale_trans_map_non_rigid, point_xi, t,
                        f_of_t_list);
    }

    t += dt;
  }

  // test serialized/deserialized map and names
  const auto rot_map_deserialized = serialize_and_deserialize(rot_map);
  const auto scale_map_deserialized =
      serialize_and_deserialize(scale_map_non_rigid);
  const auto trans_map_non_rigid_deserialized =
      serialize_and_deserialize(trans_map_non_rigid);
  const auto rot_scale_map_deserialized =
      serialize_and_deserialize(rot_scale_map_non_rigid);
  const auto rot_trans_map_non_rigid_deserialized =
      serialize_and_deserialize(rot_trans_map_non_rigid);
  const auto scale_trans_map_non_rigid_deserialized =
      serialize_and_deserialize(scale_trans_map_non_rigid);
  const auto rot_scale_trans_map_non_rigid_deserialized =
      serialize_and_deserialize(rot_scale_trans_map_non_rigid);

  CHECK(rot_map == rot_map_deserialized);
  CHECK(scale_map_non_rigid == scale_map_deserialized);
  CHECK(trans_map_non_rigid == trans_map_non_rigid_deserialized);
  CHECK(rot_scale_map_non_rigid == rot_scale_map_deserialized);
  CHECK(rot_trans_map_non_rigid == rot_trans_map_non_rigid_deserialized);
  CHECK(scale_trans_map_non_rigid == scale_trans_map_non_rigid_deserialized);
  CHECK(rot_scale_trans_map_non_rigid ==
        rot_scale_trans_map_non_rigid_deserialized);

  const auto check_names1 = [](const auto& names) {
    CHECK(names.size() == 1);
    CHECK((names.count("rotation_angle") == 1 or
           names.count("translation") == 1));
  };
  const auto check_names2 = [](const auto& names) {
    CHECK(names.size() == 2);
    CHECK((
        (names.count("rotation_angle") == 1 and
         names.count("translation") == 1) or
        (names.count("expansion_a") == 1 and names.count("expansion_b") == 1)));
  };
  const auto check_names3 = [](const auto& names) {
    CHECK(names.size() == 3);
    CHECK((
        (names.count("rotation_angle") == 1 and
         names.count("expansion_a") == 1 and names.count("expansion_b") == 1) or
        (names.count("translation") == 1 and names.count("expansion_a") == 1 and
         names.count("expansion_b") == 1)));
  };
  const auto check_names4 = [](const auto& names) {
    CHECK(names.size() == 4);
    CHECK((names.count("rotation_angle") == 1 and
           names.count("expansion_a") == 1 and
           names.count("expansion_b") == 1 and
           names.count("translation") == 1));
  };
  check_names1(rot_map.function_of_time_names());
  check_names1(trans_map_non_rigid.function_of_time_names());
  check_names2(scale_map_non_rigid.function_of_time_names());
  check_names2(rot_trans_map_non_rigid.function_of_time_names());
  check_names3(rot_scale_map_non_rigid.function_of_time_names());
  check_names3(scale_trans_map_non_rigid.function_of_time_names());
  check_names4(rot_scale_trans_map_non_rigid.function_of_time_names());
}
namespace domain {
// [[Timeout, 10]]
SPECTRE_TEST_CASE("Unit.Domain.CoordinateMaps.TimeDependent.RotScaleTrans",
                  "[Domain][Unit]") {
  test_RotScaleTrans<2>();
  test_RotScaleTrans<3>();
}
}  // namespace domain
