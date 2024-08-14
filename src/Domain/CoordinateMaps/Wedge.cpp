// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/Wedge.hpp"

#include <climits>
#include <cmath>
#include <cstddef>
#include <pup.h>

#include "DataStructures/Tensor/EagerMath/Determinant.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/Distribution.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/DereferenceWrapper.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace domain::CoordinateMaps {

template <size_t Dim>
Wedge<Dim>::Wedge(const double radius_inner, const double radius_outer,
                  const double sphericity_inner, const double sphericity_outer,
                  const double cube_half_length,
                  const std::array<double, Dim> focal_offset,
                  OrientationMap<Dim> orientation_of_wedge,
                  const bool with_equiangular_map,
                  const WedgeHalves halves_to_use,
                  const Distribution radial_distribution,
                  const std::array<double, Dim - 1>& opening_angles,
                  const bool with_adapted_equiangular_map)
    : radius_inner_(radius_inner),
      radius_outer_(radius_outer),
      sphericity_inner_(sphericity_inner),
      sphericity_outer_(sphericity_outer),
      cube_half_length_(cube_half_length),
      focal_offset_(focal_offset),
      orientation_of_wedge_(std::move(orientation_of_wedge)),
      with_equiangular_map_(with_equiangular_map),
      halves_to_use_(halves_to_use),
      radial_distribution_(radial_distribution),
      opening_angles_(opening_angles) {
  const double sqrt_dim = sqrt(double{Dim});
  const bool approx_zero_offset =
      equal_within_roundoff(magnitude(focal_offset_), 0.0);
  ASSERT(radius_inner > 0.0,
         "The radius of the inner surface must be greater than zero.");
  ASSERT(sphericity_inner >= 0.0 and sphericity_inner <= 1.0,
         "Sphericity of the inner surface must be between 0 and 1");
  ASSERT(sphericity_outer >= 0.0 and sphericity_outer <= 1.0,
         "Sphericity of the outer surface must be between 0 and 1");
  ASSERT(radius_outer > radius_inner,
         "The radius of the outer surface must be greater than the radius of "
         "the inner surface.");
  ASSERT(radial_distribution_ == Distribution::Linear or
             (sphericity_inner_ == 1.0 and sphericity_outer_ == 1.0),
         "Only the 'Linear' radial distribution is supported for non-spherical "
         "wedges.");
  ASSERT(
      get(determinant(discrete_rotation_jacobian(orientation_of_wedge_))) > 0.0,
      "Wedge rotations must be done in such a manner that the sign of "
      "the determinant of the discrete rotation is positive. This is to "
      "preserve handedness of the coordinates.");
  ASSERT(approx_zero_offset or (opening_angles_ == make_array<Dim - 1>(M_PI_2)),
         "Cannot use both a non-zero focal offset and opening angles not equal "
         "to pi/2.");
  ASSERT(opening_angles_ != make_array<Dim - 1>(M_PI_2) ? with_equiangular_map
                                                        : true,
         "If using opening angles other than pi/2, then the "
         "equiangular map option must be turned on.");

  if (approx_zero_offset) {
    ASSERT(radius_outer *
                   ((1.0 - sphericity_outer_) / sqrt_dim + sphericity_outer_) >
               radius_inner *
                   ((1.0 - sphericity_inner) / sqrt_dim + sphericity_inner),
           "The arguments passed into the constructor for Wedge result in an "
           "object where the outer surface is pierced by the inner surface.");
  } else {
    ASSERT(sphericity_inner_ == 1.0,
           "Focal offsets are not supported for inner sphericity < 1.0");
    ASSERT(sphericity_outer_ == 0.0 or sphericity_outer_ == 1.0,  // NOLINT
           "Focal offsets are only supported for wedges with outer sphericity "
           "of 1.0 or 0.0");

    // coord of focal_offset_ with largest magnitude
    const double max_abs_focal_offset_coord = *alg::max_element(
        focal_offset_,
        [](const int& a, const int& b) { return abs(a) < abs(b); });

    if (sphericity_outer_ == 1.0) {
      // note: this assert may be more restrictive than we need, can be revisted
      // if needed
      ASSERT(
          max_abs_focal_offset_coord + radius_outer_ < cube_half_length_,
          "For a spherical focally offset Wedge, the sum of the outer radius "
          "and the coordinate of the focal offset with the largest magnitude "
          "must be less than the cube half length. In other words, the "
          "spherical surface at the given outer radius centered at the focal "
          "offset must not pierce the cube of length 2 * cube_half_length_ "
          "centered at the origin. See the Wedge class documentation for a "
          "visual representation of this sphere and cube.");

    } else if (sphericity_outer_ == 0.0) {
      // if sphericity_outer_= 0.0, the outer surface of the Wedge is the parent
      // surface
      ASSERT(
          max_abs_focal_offset_coord + radius_inner_ < cube_half_length_,
          "For a cubical focally offset Wedge, the sum of the inner radius "
          "and the coordinate of the focal offset with the largest magnitude "
          "must be less than the cube half length. In other words, the "
          "spherical surface at the given inner radius centered at the focal "
          "offset must not pierce the cube of length 2 * cube_half_length_ "
          "centered at the origin. See the Wedge class documentation for a "
          "visual representation of this sphere and cube.");
    }
  }

  if (radial_distribution_ == Distribution::Linear) {
    sphere_zero_ = 0.5 * (sphericity_outer_ * radius_outer +
                          sphericity_inner * radius_inner);
    sphere_rate_ = 0.5 * (sphericity_outer_ * radius_outer -
                          sphericity_inner * radius_inner);
    if (approx_zero_offset) {
      scaled_frustum_zero_ = 0.5 / sqrt_dim *
                             ((1.0 - sphericity_outer_) * radius_outer +
                              (1.0 - sphericity_inner) * radius_inner);
      scaled_frustum_rate_ = 0.5 / sqrt_dim *
                             ((1.0 - sphericity_outer_) * radius_outer -
                              (1.0 - sphericity_inner) * radius_inner);
    } else {
      scaled_frustum_zero_ =
          0.5 * cube_half_length_ *
          ((1.0 - sphericity_outer_) + (1.0 - sphericity_inner));
      scaled_frustum_rate_ =
          0.5 * cube_half_length_ *
          ((1.0 - sphericity_outer_) - (1.0 - sphericity_inner));
    }
  } else if (radial_distribution_ == Distribution::Logarithmic) {
    scaled_frustum_zero_ = 0.0;
    sphere_zero_ = 0.5 * (log(radius_outer * radius_inner));
    scaled_frustum_rate_ = 0.0;
    sphere_rate_ = 0.5 * (log(radius_outer / radius_inner));
  } else if (radial_distribution_ == Distribution::Inverse) {
    scaled_frustum_zero_ = 0.0;
    // Most places where sphere_zero_ and sphere_rate_ would be used
    // would cause precision issues for large wedges, so we don't
    // define them.
    // 0.5 * (radius_inner + radius_outer) / radius_inner / radius_outer;
    sphere_zero_ = std::numeric_limits<double>::signaling_NaN();
    scaled_frustum_rate_ = 0.0;
    // 0.5 * (radius_inner - radius_outer) / radius_inner / radius_outer;
    sphere_rate_ = std::numeric_limits<double>::signaling_NaN();
  } else {
    ERROR("Unsupported radial distribution: " << radial_distribution_);
  }
  if (with_adapted_equiangular_map) {
    opening_angles_distribution_ = opening_angles_;
  } else {
    opening_angles_distribution_ = make_array<Dim - 1>(M_PI_2);
  }
}

template <size_t Dim>
template <typename T>
tt::remove_cvref_wrap_t<T> Wedge<Dim>::get_s_factor(const T& zeta) const {
  if (radial_distribution_ == Distribution::Linear) {
    return (sphere_zero_ + sphere_rate_ * zeta);
  } else if (radial_distribution_ == Distribution::Logarithmic) {
    return exp(sphere_zero_ + sphere_rate_ * zeta);
  } else if (radial_distribution_ == Distribution::Inverse) {
    return 2.0 / ((1.0 + zeta) / radius_outer_ + (1.0 - zeta) / radius_inner_);
  } else {
    ERROR("Unsupported radial distribution: " << radial_distribution_);
  }
}

template <size_t Dim>
template <typename T>
tt::remove_cvref_wrap_t<T> Wedge<Dim>::get_s_factor_deriv(
    const T& zeta, const T& s_factor) const {
  if (radial_distribution_ == Distribution::Linear) {
    return make_with_value<T>(zeta, sphere_rate_);
  } else if (radial_distribution_ == Distribution::Logarithmic) {
    return 0.5 * s_factor * log(radius_outer_ / radius_inner_);
  } else if (radial_distribution_ == Distribution::Inverse) {
    return 2.0 *
           ((radius_inner_ * square(radius_outer_)) -
            square(radius_inner_) * radius_outer_) /
           square(radius_inner_ + radius_outer_ +
                  zeta * (radius_inner_ - radius_outer_));
  } else {
    ERROR("Unsupported radial distribution: " << radial_distribution_);
  }
}

template <size_t Dim>
template <typename T>
tt::remove_cvref_wrap_t<T> Wedge<Dim>::get_generalized_z(
    const T& zeta, const T& one_over_rho, const T& s_factor) const {
  if (radial_distribution_ == Distribution::Linear) {
    return s_factor * one_over_rho +
           (scaled_frustum_zero_ + scaled_frustum_rate_ * zeta);
  } else if (radial_distribution_ == Distribution::Logarithmic or
             radial_distribution_ == Distribution::Inverse) {
    return s_factor * one_over_rho;
  } else {
    ERROR("Unsupported radial distribution: " << radial_distribution_);
  }
}

template <size_t Dim>
template <typename T>
tt::remove_cvref_wrap_t<T> Wedge<Dim>::get_generalized_z(
    const T& zeta, const T& one_over_rho) const {
  return get_generalized_z(zeta, one_over_rho, get_s_factor(zeta));
}

template <size_t Dim>
template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, Dim> Wedge<Dim>::get_d_generalized_z(
    const T& zeta, const T& one_over_rho, const T& s_factor,
    const std::array<tt::remove_cvref_wrap_t<T>, Dim - 1>& cap_deriv,
    const std::array<tt::remove_cvref_wrap_t<T>, Dim>& rho_vec) const {
  using ReturnType = tt::remove_cvref_wrap_t<T>;

  const ReturnType one_over_rho_cubed = pow<3>(one_over_rho);
  const ReturnType s_factor_over_rho_cubed = s_factor * one_over_rho_cubed;

  std::array<ReturnType, Dim> d_generalized_z{};
  // Polar angle
  d_generalized_z[polar_coord] =
      -s_factor_over_rho_cubed * cap_deriv[0] * rho_vec[polar_coord];
  // Radial coordinate
  if (radial_distribution_ == Distribution::Linear) {
    // note: sphere_rate_ = s_factor_deriv for Linear
    d_generalized_z[radial_coord] =
        sphere_rate_ * one_over_rho + scaled_frustum_rate_;
  } else if (radial_distribution_ == Distribution::Logarithmic or
             radial_distribution_ == Distribution::Inverse) {
    const ReturnType s_factor_deriv = get_s_factor_deriv(zeta, s_factor);
    d_generalized_z[radial_coord] = s_factor_deriv * one_over_rho;
  } else {
    ERROR("Unsupported radial distribution: " << radial_distribution_);
  }
  if constexpr (Dim == 3) {
    // Azimuthal angle
    d_generalized_z[azimuth_coord] =
        -s_factor_over_rho_cubed * cap_deriv[1] * rho_vec[azimuth_coord];
  }

  return d_generalized_z;
}

template <size_t Dim>
template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, Dim> Wedge<Dim>::operator()(
    const std::array<T, Dim>& source_coords) const {
  using ReturnType = tt::remove_cvref_wrap_t<T>;

  // Radial coordinate
  const ReturnType& zeta = source_coords[radial_coord];

  // Polar angle
  ReturnType xi = source_coords[polar_coord];
  if (halves_to_use_ == WedgeHalves::UpperOnly) {
    xi += 1.0;
    xi *= 0.5;
  } else if (halves_to_use_ == WedgeHalves::LowerOnly) {
    xi -= 1.0;
    xi *= 0.5;
  }

  const auto rotated_focus =
      discrete_rotation(orientation_of_wedge_.inverse_map(), focal_offset_);

  std::array<ReturnType, Dim - 1> cap{};
  cap[0] = with_equiangular_map_
               ? tan(0.5 * opening_angles_[0]) *
                     tan(0.5 * opening_angles_distribution_[0] * xi) /
                     tan(0.5 * opening_angles_distribution_[0])
               : xi;
  ReturnType one_over_rho =
      square(1.0 - rotated_focus[radial_coord] / cube_half_length_) +
      square(cap[0] - rotated_focus[polar_coord] / cube_half_length_);
  if constexpr (Dim == 3) {
    // Azimuthal angle
    const ReturnType& eta = source_coords[azimuth_coord];
    cap[1] = with_equiangular_map_
                 ? tan(0.5 * opening_angles_[1]) *
                       tan(0.5 * opening_angles_distribution_[1] * eta) /
                       tan(0.5 * opening_angles_distribution_[1])
                 : eta;
    one_over_rho +=
        square(cap[1] - rotated_focus[azimuth_coord] / cube_half_length_);
  }
  one_over_rho = 1.0 / sqrt(one_over_rho);

  const ReturnType generalized_z = get_generalized_z(zeta, one_over_rho);
  std::array<ReturnType, Dim> physical_coords{};
  physical_coords[radial_coord] =
      generalized_z * (1.0 - rotated_focus[radial_coord] / cube_half_length_) +
      rotated_focus[radial_coord];
  physical_coords[polar_coord] =
      generalized_z *
          (cap[0] - rotated_focus[polar_coord] / cube_half_length_) +
      rotated_focus[polar_coord];
  if constexpr (Dim == 3) {
    physical_coords[azimuth_coord] =
        generalized_z *
            (cap[1] - rotated_focus[azimuth_coord] / cube_half_length_) +
        rotated_focus[azimuth_coord];
  }
  return discrete_rotation(orientation_of_wedge_, std::move(physical_coords));
}

template <size_t Dim>
std::optional<std::array<double, Dim>> Wedge<Dim>::inverse(
    const std::array<double, Dim>& target_coords) const {
  const std::array<double, Dim> physical_coords =
      discrete_rotation(orientation_of_wedge_.inverse_map(), target_coords);
  const auto rotated_focus =
      discrete_rotation(orientation_of_wedge_.inverse_map(), focal_offset_);

  if (physical_coords[radial_coord] < rotated_focus[radial_coord] or
      equal_within_roundoff(physical_coords[radial_coord],
                            rotated_focus[radial_coord])) {
    return std::nullopt;
  }

  const double generalized_z =
      (physical_coords[radial_coord] - rotated_focus[radial_coord]) /
      (1.0 - rotated_focus[radial_coord] / cube_half_length_);
  const double one_over_generalized_z = 1.0 / generalized_z;

  std::array<double, Dim - 1> cap{};
  cap[0] = (physical_coords[polar_coord] - rotated_focus[polar_coord]) *
               one_over_generalized_z +
           rotated_focus[polar_coord] / cube_half_length_;
  if constexpr (Dim == 3) {
    cap[1] = (physical_coords[azimuth_coord] - rotated_focus[azimuth_coord]) *
                 one_over_generalized_z +
             rotated_focus[azimuth_coord] / cube_half_length_;
  }

  // Radial coordinate
  double zeta = std::numeric_limits<double>::signaling_NaN();
  const double radius = magnitude(physical_coords - rotated_focus);
  if (radial_distribution_ == Distribution::Linear) {
    const double one_over_rho = generalized_z / radius;
    const double zeta_coefficient =
        (scaled_frustum_rate_ + sphere_rate_ * one_over_rho);
    // If -sphere_rate_/scaled_frustum_rate_ > 1, then
    // there exists a cone in x,y,z space given by the surface
    // zeta_coefficient=0; the map is singular on this surface.
    // We return nullopt if we are on or outside this cone.
    // If scaled_frustum_rate_ > 0, then outside the cone
    // corresponds to zeta_coefficient > 0, and if scaled_frustum_rate_
    // < 0, then outside the cone corresponds to zeta_coefficient < 0.
    // We test in two cases, and avoid division.
    if ((scaled_frustum_rate_ > 0.0 and scaled_frustum_rate_ < -sphere_rate_ and
         zeta_coefficient > 0.0) or
        (scaled_frustum_rate_ < 0.0 and scaled_frustum_rate_ > -sphere_rate_ and
         zeta_coefficient < 0.0) or
        equal_within_roundoff(zeta_coefficient, 0.0)) {
      return std::nullopt;
    }
    const auto z_zero = (scaled_frustum_zero_ + sphere_zero_ * one_over_rho);
    zeta = (generalized_z - z_zero) / zeta_coefficient;
  } else if (radial_distribution_ == Distribution::Logarithmic) {
    zeta = (log(radius) - sphere_zero_) / sphere_rate_;
  } else if (radial_distribution_ == Distribution::Inverse) {
    zeta = (radius_inner_ * (radius_outer_ / radius - 1.0) +
            radius_outer_ * (radius_inner_ / radius - 1.0)) /
           (radius_inner_ - radius_outer_);
  } else {
    ERROR("Unsupported radial distribution: " << radial_distribution_);
  }
  // Polar angle
  double xi = with_equiangular_map_
                  ? 2.0 *
                        atan(tan(0.5 * opening_angles_distribution_[0]) /
                             tan(0.5 * opening_angles_[0]) * cap[0]) /
                        opening_angles_distribution_[0]
                  : cap[0];
  if (halves_to_use_ == WedgeHalves::UpperOnly) {
    xi *= 2.0;
    xi -= 1.0;
  } else if (halves_to_use_ == WedgeHalves::LowerOnly) {
    xi *= 2.0;
    xi += 1.0;
  }
  std::array<double, Dim> logical_coords{};
  logical_coords[radial_coord] = zeta;
  logical_coords[polar_coord] = xi;
  if constexpr (Dim == 3) {
    // Azimuthal angle
    logical_coords[azimuth_coord] =
        with_equiangular_map_
            ? 2.0 *
                  atan(tan(0.5 * opening_angles_distribution_[1]) /
                       tan(0.5 * opening_angles_[1]) * cap[1]) /
                  opening_angles_distribution_[1]
            : cap[1];
  }
  return logical_coords;
}

template <size_t Dim>
template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame> Wedge<Dim>::jacobian(
    const std::array<T, Dim>& source_coords) const {
  using ReturnType = tt::remove_cvref_wrap_t<T>;

  // Radial coordinate
  const ReturnType& zeta = source_coords[radial_coord];

  // Polar angle
  ReturnType xi = source_coords[polar_coord];
  if (halves_to_use_ == WedgeHalves::UpperOnly) {
    xi += 1.0;
    xi *= 0.5;
  } else if (halves_to_use_ == WedgeHalves::LowerOnly) {
    xi -= 1.0;
    xi *= 0.5;
  }

  std::array<ReturnType, Dim - 1> cap{};
  std::array<ReturnType, Dim - 1> cap_deriv{};
  cap[0] = with_equiangular_map_
               ? tan(0.5 * opening_angles_[0]) *
                     tan(0.5 * opening_angles_distribution_[0] * xi) /
                     tan(0.5 * opening_angles_distribution_[0])
               : xi;
  cap_deriv[0] =
      with_equiangular_map_
          ? 0.5 * opening_angles_distribution_[0] *
                tan(0.5 * opening_angles_[0]) /
                tan(0.5 * opening_angles_distribution_[0]) /
                square(cos(0.5 * opening_angles_distribution_[0] * xi))
          : make_with_value<ReturnType>(xi, 1.0);

  const auto rotated_focus =
      discrete_rotation(orientation_of_wedge_.inverse_map(), focal_offset_);

  std::array<ReturnType, Dim> rho_vec{};
  rho_vec[polar_coord] =
      cap[0] - rotated_focus[polar_coord] / cube_half_length_;
  rho_vec[radial_coord] = make_with_value<ReturnType>(cap[0], 1.0) -
                          rotated_focus[radial_coord] / cube_half_length_;

  ReturnType one_over_rho =
      square(1.0 - rotated_focus[radial_coord] / cube_half_length_) +
      square(cap[0] - rotated_focus[polar_coord] / cube_half_length_);

  if constexpr (Dim == 3) {
    // Azimuthal angle
    const ReturnType& eta = source_coords[azimuth_coord];
    cap[1] = with_equiangular_map_
                 ? tan(0.5 * opening_angles_[1]) *
                       tan(0.5 * opening_angles_distribution_[1] * eta) /
                       tan(0.5 * opening_angles_distribution_[1])
                 : eta;
    cap_deriv[1] =
        with_equiangular_map_
            ? 0.5 * opening_angles_distribution_[1] *
                  tan(0.5 * opening_angles_[1]) /
                  tan(0.5 * opening_angles_distribution_[1]) /
                  square(cos(0.5 * opening_angles_distribution_[1] * eta))
            : make_with_value<ReturnType>(xi, 1.0);

    rho_vec[azimuth_coord] =
        cap[1] - rotated_focus[azimuth_coord] / cube_half_length_;

    one_over_rho +=
        square(cap[1] - rotated_focus[azimuth_coord] / cube_half_length_);
  }
  one_over_rho = 1.0 / sqrt(one_over_rho);

  const ReturnType s_factor = get_s_factor(zeta);
  const ReturnType generalized_z =
      get_generalized_z(zeta, one_over_rho, s_factor);
  const std::array<ReturnType, Dim> d_generalized_z =
      get_d_generalized_z(zeta, one_over_rho, s_factor, cap_deriv, rho_vec);

  auto jacobian_matrix =
      make_with_value<tnsr::Ij<ReturnType, Dim, Frame::NoFrame>>(xi, 0.0);

  // Derivative by polar angle
  std::array<ReturnType, Dim> dxyz_dxi{};
  dxyz_dxi[radial_coord] = rho_vec[radial_coord] * d_generalized_z[polar_coord];
  dxyz_dxi[polar_coord] = rho_vec[polar_coord] * d_generalized_z[polar_coord] +
                          cap_deriv[0] * generalized_z;

  if constexpr (Dim == 3) {
    dxyz_dxi[azimuth_coord] =
        rho_vec[azimuth_coord] * d_generalized_z[polar_coord];
  }

  // Implement Scalings:
  if (halves_to_use_ != WedgeHalves::Both) {
    for (size_t d = 0; d < Dim; ++d) {
      gsl::at(dxyz_dxi, d) *= 0.5;
    }
  }

  std::array<ReturnType, Dim> dX_dlogical =
      discrete_rotation(orientation_of_wedge_, std::move(dxyz_dxi));
  get<0, polar_coord>(jacobian_matrix) = dX_dlogical[0];
  get<1, polar_coord>(jacobian_matrix) = dX_dlogical[1];
  if constexpr (Dim == 3) {
    get<2, polar_coord>(jacobian_matrix) = dX_dlogical[2];
  }

  // Derivative by azimuthal angle
  if constexpr (Dim == 3) {
    std::array<ReturnType, Dim> dxyz_deta{};
    dxyz_deta[radial_coord] =
        rho_vec[radial_coord] * d_generalized_z[azimuth_coord];

    dxyz_deta[azimuth_coord] =
        rho_vec[azimuth_coord] * d_generalized_z[azimuth_coord] +
        cap_deriv[1] * generalized_z;

    dxyz_deta[polar_coord] =
        rho_vec[polar_coord] * d_generalized_z[azimuth_coord];

    dX_dlogical =
        discrete_rotation(orientation_of_wedge_, std::move(dxyz_deta));
    get<0, azimuth_coord>(jacobian_matrix) = dX_dlogical[0];
    get<1, azimuth_coord>(jacobian_matrix) = dX_dlogical[1];
    get<2, azimuth_coord>(jacobian_matrix) = dX_dlogical[2];
  }

  // Derivative by radial coordinate
  std::array<ReturnType, Dim> dxyz_dzeta{};
  dxyz_dzeta[radial_coord] =
      rho_vec[radial_coord] * d_generalized_z[radial_coord];
  dxyz_dzeta[polar_coord] =
      rho_vec[polar_coord] * d_generalized_z[radial_coord];

  if constexpr (Dim == 3) {
    dxyz_dzeta[azimuth_coord] =
        rho_vec[azimuth_coord] * d_generalized_z[radial_coord];
  }

  dX_dlogical = discrete_rotation(orientation_of_wedge_, std::move(dxyz_dzeta));
  get<0, radial_coord>(jacobian_matrix) = dX_dlogical[0];
  get<1, radial_coord>(jacobian_matrix) = dX_dlogical[1];
  if constexpr (Dim == 3) {
    get<2, radial_coord>(jacobian_matrix) = dX_dlogical[2];
  }

  return jacobian_matrix;
}

template <size_t Dim>
template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame>
Wedge<Dim>::inv_jacobian(const std::array<T, Dim>& source_coords) const {
  using ReturnType = tt::remove_cvref_wrap_t<T>;

  // Radial coordinate
  const ReturnType& zeta = source_coords[radial_coord];

  // Polar angle
  ReturnType xi = source_coords[polar_coord];
  if (halves_to_use_ == WedgeHalves::UpperOnly) {
    xi += 1.0;
    xi *= 0.5;
  } else if (halves_to_use_ == WedgeHalves::LowerOnly) {
    xi -= 1.0;
    xi *= 0.5;
  }

  std::array<ReturnType, Dim - 1> cap{};
  std::array<ReturnType, Dim - 1> cap_deriv{};
  cap[0] = with_equiangular_map_
               ? tan(0.5 * opening_angles_[0]) *
                     tan(0.5 * opening_angles_distribution_[0] * xi) /
                     tan(0.5 * opening_angles_distribution_[0])
               : xi;
  cap_deriv[0] =
      with_equiangular_map_
          ? 0.5 * opening_angles_distribution_[0] *
                tan(0.5 * opening_angles_[0]) /
                tan(0.5 * opening_angles_distribution_[0]) /
                square(cos(0.5 * opening_angles_distribution_[0] * xi))
          : make_with_value<ReturnType>(xi, 1.0);

  const auto rotated_focus =
      discrete_rotation(orientation_of_wedge_.inverse_map(), focal_offset_);

  std::array<ReturnType, Dim> rho_vec{};
  rho_vec[polar_coord] =
      cap[0] - rotated_focus[polar_coord] / cube_half_length_;
  rho_vec[radial_coord] = make_with_value<ReturnType>(cap[0], 1.0) -
                          rotated_focus[radial_coord] / cube_half_length_;

  ReturnType one_over_rho =
      square(1.0 - rotated_focus[radial_coord] / cube_half_length_) +
      square(cap[0] - rotated_focus[polar_coord] / cube_half_length_);

  if constexpr (Dim == 3) {
    // Azimuthal angle
    const ReturnType& eta = source_coords[azimuth_coord];
    cap[1] = with_equiangular_map_
                 ? tan(0.5 * opening_angles_[1]) *
                       tan(0.5 * opening_angles_distribution_[1] * eta) /
                       tan(0.5 * opening_angles_distribution_[1])
                 : eta;
    cap_deriv[1] =
        with_equiangular_map_
            ? 0.5 * opening_angles_distribution_[1] *
                  tan(0.5 * opening_angles_[1]) /
                  tan(0.5 * opening_angles_distribution_[1]) /
                  square(cos(0.5 * opening_angles_distribution_[1] * eta))
            : make_with_value<ReturnType>(xi, 1.0);

    rho_vec[azimuth_coord] =
        cap[1] - rotated_focus[azimuth_coord] / cube_half_length_;

    one_over_rho +=
        square(cap[1] - rotated_focus[azimuth_coord] / cube_half_length_);
  }
  one_over_rho = 1.0 / sqrt(one_over_rho);

  const ReturnType s_factor = get_s_factor(zeta);
  const ReturnType generalized_z =
      get_generalized_z(zeta, one_over_rho, s_factor);
  const ReturnType one_over_generalized_z = 1.0 / generalized_z;
  const std::array<ReturnType, Dim> d_generalized_z =
      get_d_generalized_z(zeta, one_over_rho, s_factor, cap_deriv, rho_vec);
  const ReturnType one_over_d_generalized_z_dzeta =
      1.0 / d_generalized_z[radial_coord];
  const ReturnType one_over_rho_z = 1.0 / rho_vec[radial_coord];
  const ReturnType scaled_z_frustum =
      scaled_frustum_zero_ + scaled_frustum_rate_ * zeta;

  auto inv_jacobian_matrix =
      make_with_value<tnsr::Ij<ReturnType, Dim, Frame::NoFrame>>(xi, 0.0);

  // Derivatives of polar angle
  std::array<ReturnType, Dim> dxi_dxyz{};
  dxi_dxyz[polar_coord] = 1.0 / (generalized_z * cap_deriv[0]);
  // Implement Scalings:
  if (halves_to_use_ != WedgeHalves::Both) {
    dxi_dxyz[polar_coord] *= 2.0;
  }

  dxi_dxyz[radial_coord] =
      -dxi_dxyz[polar_coord] * one_over_rho_z * rho_vec[polar_coord];

  if constexpr (Dim == 3) {
    dxi_dxyz[azimuth_coord] = make_with_value<ReturnType>(xi, 0.0);
  }

  std::array<ReturnType, Dim> dlogical_dX =
      discrete_rotation(orientation_of_wedge_, std::move(dxi_dxyz));
  get<polar_coord, 0>(inv_jacobian_matrix) = dlogical_dX[0];
  get<polar_coord, 1>(inv_jacobian_matrix) = dlogical_dX[1];
  if constexpr (Dim == 3) {
    get<polar_coord, 2>(inv_jacobian_matrix) = dlogical_dX[2];
  }

  // Derivatives of radial coordinate

  // a common term that appears in the Jacobian, see Wedge docs
  const ReturnType T_factor =
      s_factor * one_over_d_generalized_z_dzeta * pow<3>(one_over_rho);

  std::array<ReturnType, Dim> dzeta_dxyz{};
  dzeta_dxyz[polar_coord] =
      T_factor * rho_vec[polar_coord] * one_over_generalized_z;
  dzeta_dxyz[radial_coord] =
      one_over_generalized_z *
      (one_over_rho_z * scaled_z_frustum * one_over_d_generalized_z_dzeta +
       T_factor * rho_vec[radial_coord]);
  if constexpr (Dim == 3) {
    dzeta_dxyz[azimuth_coord] =
        T_factor * rho_vec[azimuth_coord] * one_over_generalized_z;
  }

  dlogical_dX = discrete_rotation(orientation_of_wedge_, std::move(dzeta_dxyz));
  get<radial_coord, 0>(inv_jacobian_matrix) = dlogical_dX[0];
  get<radial_coord, 1>(inv_jacobian_matrix) = dlogical_dX[1];
  if constexpr (Dim == 3) {
    get<radial_coord, 2>(inv_jacobian_matrix) = dlogical_dX[2];
  }

  if constexpr (Dim == 3) {
    // Derivatives of azimuthal angle
    std::array<ReturnType, Dim> deta_dxyz{};
    deta_dxyz[polar_coord] = make_with_value<ReturnType>(xi, 0.0);
    deta_dxyz[azimuth_coord] = 1.0 / (generalized_z * cap_deriv[1]);
    deta_dxyz[radial_coord] =
        -deta_dxyz[azimuth_coord] * one_over_rho_z * rho_vec[azimuth_coord];

    dlogical_dX =
        discrete_rotation(orientation_of_wedge_, std::move(deta_dxyz));
    get<azimuth_coord, 0>(inv_jacobian_matrix) = dlogical_dX[0];
    get<azimuth_coord, 1>(inv_jacobian_matrix) = dlogical_dX[1];
    get<azimuth_coord, 2>(inv_jacobian_matrix) = dlogical_dX[2];
  }

  return inv_jacobian_matrix;
}

template <size_t Dim>
void Wedge<Dim>::pup(PUP::er& p) {
  size_t version = 2;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | radius_inner_;
    p | radius_outer_;
    p | sphericity_inner_;
    p | sphericity_outer_;
    p | orientation_of_wedge_;
    p | with_equiangular_map_;
    p | halves_to_use_;
    p | radial_distribution_;
    p | scaled_frustum_zero_;
    p | sphere_zero_;
    p | scaled_frustum_rate_;
    p | sphere_rate_;
    if (version == 0) {
      double half_opening_angle = std::numeric_limits<double>::signaling_NaN();
      p | half_opening_angle;
      opening_angles_[0] = 2. * half_opening_angle;
      if constexpr (Dim == 3) {
        opening_angles_[1] = M_PI_2;
      }
      opening_angles_distribution_ = opening_angles_;
    }
  }
  if (version >= 1) {
    p | opening_angles_;
    p | opening_angles_distribution_;
  }
  if (version >= 2) {
    p | cube_half_length_;
    p | focal_offset_;
  } else if (p.isUnpacking()) {
    // While the length of the cube outside of an old Wedge is not known
    // internally to the Wedge class, any length will mathematically work
    // because all old Wedges have no focal offset and any terms in the math
    // that use the cube length will be zeroed out by the zero offset. We choose
    // this half length to be 1.0 just because it's a simple value.
    cube_half_length_ = 1.0;
    focal_offset_ = make_with_value<std::array<double, Dim>>(Dim, 0.0);
  }
}

template <size_t Dim>
bool operator==(const Wedge<Dim>& lhs, const Wedge<Dim>& rhs) {
  return lhs.radius_inner_ == rhs.radius_inner_ and
         lhs.radius_outer_ == rhs.radius_outer_ and
         lhs.cube_half_length_ == rhs.cube_half_length_ and
         lhs.focal_offset_ == rhs.focal_offset_ and
         lhs.orientation_of_wedge_ == rhs.orientation_of_wedge_ and
         lhs.with_equiangular_map_ == rhs.with_equiangular_map_ and
         lhs.halves_to_use_ == rhs.halves_to_use_ and
         lhs.radial_distribution_ == rhs.radial_distribution_ and
         lhs.sphericity_inner_ == rhs.sphericity_inner_ and
         lhs.sphericity_outer_ == rhs.sphericity_outer_ and
         lhs.scaled_frustum_zero_ == rhs.scaled_frustum_zero_ and
         lhs.scaled_frustum_rate_ == rhs.scaled_frustum_rate_ and
         lhs.opening_angles_ == rhs.opening_angles_ and
         lhs.opening_angles_distribution_ == rhs.opening_angles_distribution_;
}

template <size_t Dim>
bool operator!=(const Wedge<Dim>& lhs, const Wedge<Dim>& rhs) {
  return not(lhs == rhs);
}

// Explicit instantiations
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_DIM(_, data)                         \
  template class Wedge<DIM(data)>;                       \
  template bool operator==(const Wedge<DIM(data)>& lhs,  \
                           const Wedge<DIM(data)>& rhs); \
  template bool operator!=(const Wedge<DIM(data)>& lhs,  \
                           const Wedge<DIM(data)>& rhs);

#define INSTANTIATE_DTYPE(_, data)                                     \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)> \
  Wedge<DIM(data)>::operator()(                                        \
      const std::array<DTYPE(data), DIM(data)>& source_coords) const;  \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),   \
                    Frame::NoFrame>                                    \
  Wedge<DIM(data)>::jacobian(                                          \
      const std::array<DTYPE(data), DIM(data)>& source_coords) const;  \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),   \
                    Frame::NoFrame>                                    \
  Wedge<DIM(data)>::inv_jacobian(                                      \
      const std::array<DTYPE(data), DIM(data)>& source_coords) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_DIM, (2, 3))
GENERATE_INSTANTIATIONS(INSTANTIATE_DTYPE, (2, 3),
                        (double, DataVector,
                         std::reference_wrapper<const double>,
                         std::reference_wrapper<const DataVector>))

#undef DIM
#undef DTYPE
#undef INSTANTIATE_DIM
#undef INSTANTIATE_DTYPE
}  // namespace domain::CoordinateMaps
