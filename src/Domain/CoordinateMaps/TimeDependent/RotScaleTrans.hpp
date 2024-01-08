// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/MathFunctions/MathFunction.hpp"
#include "Utilities/TypeTraits/RemoveReferenceWrapper.hpp"

/// \cond
namespace domain::FunctionsOfTime {
class FunctionOfTime;
}  // namespace domain::FunctionsOfTime
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace domain::CoordinateMaps::TimeDependent {
/*!
 * \ingroup CoordMapsTimeDependentGroup
 * \brief RotScaleTrans map which applies a combination of rotation, expansion,
 * and translation based on which constructor is used.
 *
 * \details This map adds a rotation, expansion and translation based on what
 * types of maps are needed. Translation and expansion have piecewise functions
 * that map the coordinates $\vec{\xi}$ based on what region
 * $|\vec{\xi}|$ is in. Coordinates within the inner radius are translated by
 * the translation function of time T(t) and expanded by the inner expansion
 * function of time E_{a}(t). Coordinates in between the inner radius and outer
 * radius have a linear radial falloff applied to them. Coordinates at or beyond
 * the outer radius have no translation applied and are expanded by the outer
 * expansion function of time E_{b}(t). This map assumes that the center of the
 * map is at (0., 0., 0.). If translation or expansion is specified, then an
 * inner and outer radius must be specified as well. There is a rigid flag
 * which when set to true treats the coordinates as if they're within the inner
 * radius. This can be used as an optimization if you know which blocks in the
 * domain will always be within the inner radius.
 *
 * ## Mapped Coordinates
 * The RotScaleTrans map takes the coordinates $\vec{\xi}$ to the target
 * coordinates $\vec{\bar{\xi}}$ through
 * \f{equation}{
 * \vec{\bar{\xi}} = \left\{\begin{array}{ll}R(t)E_{a}(t)\vec{\xi} + \vec{T}(t),
 * & |\vec{\xi}| \leq R_{in}, \\ R(t)(w_{E} + E_{a}(t))\vec{\xi} + T(t) +
 * w_{T}T(t), & R_{in} < |\vec{\xi}| < 0.9 R_{out}, \\ R(t)(w_{E} +
 * E_{b}(t))\vec{\xi} + w_{T}T(t), & 0.9 R_{out} < |\vec{\xi}| < R_{out},
 * \\ R(t)E_{b}(t)\vec{\xi}, &
 * |\vec{\xi}| \geq R_{out} \end{array}\right.
 * \f}
 *
 * Where $R_{in}$ is the inner radius, $R_{out}$ is the outer radius, and
 * $w_{T}$ is the translation falloff factor and $w_{E}$ is the expansion
 * falloff factor found through
 * \f{equation}{
 * w_{E} = \left\{\begin{array}{ll}\frac{R_{out}(R_{in} - |\vec{\xi}|)(E_{a}(t)
 * - E_{b}(t))}{|\vec{\xi}|(R_{out} - R_{in})}, & R_{in} < |\vec{\xi}| < 0.9
 * R_{out}, \\ \frac{R_{in}(R_{out} - |\vec{\xi}|)(E_{a}(t) -
 * E_{b}(t))}{|\vec{\xi}|(R_{out} - R_{in})}, & 0.9 R_{out} \leq |\vec{\xi}| <
 * R_{out} \end{array}\right.
 * \f}
 *
 * and
 *
 * \f{equation}{
 * w_{T} = \left\{\begin{array}{ll}\frac{R_{in} - |\vec{\xi}|}{R_{out} -
 * R_{in}}, & R_{in} < |\vec{\xi}| < 0.9 R_{out}, \\ \frac{R_{out} -
 * |\vec{\xi}|}{R_{out} - R_{in}}, & 0.9 R_{out} \leq |\vec{\xi}| < R_{out}
 * \end{array}\right.
 * \f}
 *
 * $w_{E}$ and $w_{T}$ are calculated differently based on if you're closer
 * to the inner radius or outer radius to reduce roundoff error.
 *
 * ## Inverse
 * The inverse function maps the coordinates $\vec{\bar{\xi}}$ back to the
 * original coordinates $\vec{\xi}$ through different equations based on
 * which maps are supplied.
 *
 * ### If Rotation, Expansion, and Translation Maps are supplied then the
 * inverse is given by
 *
 * \f{equation}{
 * \vec{\xi} = \left\{\begin{array}{ll}R^{T}(t)(\frac{(\vec{\bar{\xi}} -
 * \vec{T}(t))}{E_{a}(t)}), & |\vec{\xi}| \leq R_{in}E_{a}(t),
 * \\ R^{T}(t)\frac{\vec{\bar{\xi}} - w_{T}T(t)}{w_{E}}, & R_{in} < |\vec{\xi}|
 * < 0.5 R_{out}, \\ R^{T}(t)\frac{\vec{\bar{\xi}} - (1.0 - w_{T})T(t)}{w_{E}},
 * & 0.5 R_{out} \leq |\vec{\xi}| < R_{out},
 * \\ R^{T}(t)\frac{\vec{\bar{\xi}}}{E_{b}(t)}, & |\vec{\xi}| \geq R_{out}
 * \end{array}\right.
 * \f}
 *
 * Where $w_{T}$ and $w_{E}$ are found through different quadratic solves.
 *
 * When closer to $R_{in}$ the quadratic has the form
 * \f{equation}{
 * w_{T}^2((E_{a}(t)R_{in} - E_{b}(t)R_{out})^2 - T(t)^2) +
 * 2w_{T}(E_{b}(t)R_{out}(E_{a}(t)R_{in} - E_{b}(t)R_{out}) +
 * T(t)\vec{\bar{\xi}})
 * + E_{b}(t)^2 R_{out}^2 - \vec{\bar{\xi}}^2
 * \f}
 *
 * where $w_{E} = \frac{w_{T}R_{in}E_{a}(t) + (1.0 -
 * w_{T})R_{out}E_{b}(t)}{w_{T}R_{in} + R_{out}(1.0 - w_{T})}$
 *
 * When closer to $R_{out}$ the quadratic has the form
 * \f{equation}{
 * w_{T}^2((E_{a}(t)R_{in} - E_{b}(t)R_{out})^2 - T(t)^2) +
 * 2w_{T}(E_{a}(t)R_{in}(E_{b}(t)R_{out} - E_{a}(t)R_{in}) + T(t)(T(t) -
 * \vec{\bar{\xi}})) + T(t)\vec{\bar{\xi}})
 * + E_{a}(t)^2 R_{in}^2 - (\vec{\bar{\xi}} - T(t))^2
 * \f}
 *
 * where $w_{E} = \frac{(1.0 - w_{T})R_{in}E_{a}(t) +
 * R_{out}E_{b}(t)}{R_{in}(1.0 - w_{T}) + R_{out}}$
 *
 * ### If Rotation and Expansion are supplied then the inverse is given by
 *
 * \f{equation}{
 * \vec{\xi} =
 * \left\{\begin{array}{ll}R^{T}(t)(\frac{\vec{\bar{\xi}}}{E_{a}(t)}), &
 * |\vec{\xi}| \leq R_{in}E_{a}(t),
 * \\ R^{T}(t)\frac{\vec{\bar{\xi}}}{w_{E}}, & R_{in} < |\vec{\xi}|
 * < 0.5 R_{out}, \\ R^{T}(t)\frac{\vec{\bar{\xi}}}{w_{E}},
 * & 0.5 R_{out} \leq |\vec{\xi}| < R_{out},
 * \\ R^{T}(t)\frac{\vec{\bar{\xi}}}{E_{b}(t)}, & |\vec{\xi}| \geq R_{out}
 * \end{array}\right.
 * \f}
 *
 * Where $w_{E}$ is found through different quadratic solves.
 *
 * When closer to $R_{in}$ the quadratic has the form
 * \f{equation}{
 * w^2(\vec{\bar{\xi}}^2(R_{out} - R_{in})^2 - \rho^2 R_{out}^2) +
 * 2w\rho(\bar{\xi}^2(R_{out} - R_{in}) - \rho R_{out}^2 E_{b}(t))
 * + \rho^2(\bar{\xi}^2 - (E_{b}(t) R_{out})^2)
 * \f}
 * with $\rho = R_{in}(E_{a}(t) - E_{b}(t))$
 *
 * When closer to $R_{out}$ the quadratic has the form
 * \f{equation}{
 * w^2(\vec{\bar{\xi}}^2(R_{out} - R_{in})^2 - \rho^2 R_{in}^2) +
 * 2w\rho(\bar{\xi}^2(R_{in} - R_{out}) - \rho R_{in}^2 E_{a}(t))
 * + \rho^2(\bar{\xi}^2 - (E_{a}(t) R_{in})^2)
 * \f}
 * with $\rho = R_{out}(E_{a}(t) - E_{b}(t))$
 *
 * ### If Rotation and Translation are supplied, then the inverse is given by
 *
 * \f{equation}{
 * \vec{\xi} = \left\{\begin{array}{ll}R^{T}(t)(\vec{\bar{\xi}} -
 * \vec{T}(t)), & |\vec{\xi}| \leq R_{in}E_{a}(t),
 * \\ R^{T}(t)(\vec{\bar{\xi}} - w_{T}T(t)), & R_{in} < |\vec{\xi}|
 * < 0.5 R_{out}, \\ R^{T}(t)(\vec{\bar{\xi}} - (1.0 - w_{T})T(t)),
 * & 0.5 R_{out} \leq |\vec{\xi}| < R_{out},
 * \\ R^{T}(t)\vec{\bar{\xi}}, & |\vec{\xi}| \geq R_{out}
 * \end{array}\right.
 * \f}
 *
 * Where $w_{T}$ is found through different quadratic solves.
 *
 * When closer to $R_{in}$ the quadratic has the form
 * \f{equation}{
 * w^2(T(t)^2 - (R_{out} - R_{in})^2) +
 * 2w(R_{out}(R_{out} - R_{in}) - T(t)\vec{\bar{\xi}})
 * + \vec{\bar{\xi}}^2 - R_{out}^2
 * \f}
 *
 * When closer to $R_{out}$ the quadratic has the form
 * \f{equation}{
 * w^2(T(t)^2 - (R_{out} - R_{in})^2) -
 * 2w(T(t)(T(t) - \vec{\bar{\xi}}) + R_{in}(R_{out} - R_{in})
 * + (T(t) - \vec{\bar{\xi}})^2 - R_{in}^2
 * \f}
 *
 * \note For all the maps with rotation, the inverse of rotation is the
 * transpose of the original rotation. For maps with translation, the inverse
 * map also assumes that if $\vec{\bar{\xi}} - \vec{T}(t) \leq R_{in}$ then the
 * translated point originally came from within the inner radius so it'll be
 * translated back without a quadratic solve.
 *
 * ## Frame Velocity
 * The Frame Velocity is found through different equations based on which maps
 * are supplied.
 *
 * ### If Rotation, Expansion, and Translation are supplied then the frame
 * velocity is found through
 *
 * \f{equation}{
 * \vec{v} = \left\{\begin{array}{ll}\vec{\xi}(dR(t)E_{a}(t) +
 * R(t)dE_{a}(t)) + \vec{dT}(t), & |\vec{\xi}| \leq R_{in},
 * \\ \vec{\xi}(dR(t)(E_{a}(t) + w_{E}) + R(t)(dE_{a}(t) + dw_{E})) + T(t) +
 * w_{T}T(t), & R_{in} < |\vec{\xi}| < 0.9 R_{out},
 * \\ \vec{\xi}(dR(t)(E_{b}(t) + w_{E}) + R(t)(dE_{b}(t) + dw_{E})) + w_{T}T(t),
 * & 0.9 R_{out} \leq |\vec{\xi}| < R_{out}, \\ \vec{\xi}(dR(t)E_{b}(t) +
 * R(t)dE_{b}(t)), & |\vec{\xi}| \geq R_{out} \end{array}\right.
 * \f}
 *
 * where $dw_{E}$ is the derivative of the $w_{E}$ given by
 *
 * \f{equation}{
 * dw_{E} = \left\{\begin{array}{ll}\frac{R_{out}(R_{in} -
 * |\vec{\xi}|)(dE_{a}(t)
 * - dE_{b}(t))}{|\vec{\xi}|(R_{out} - R_{in})}, & R_{in} < |\vec{\xi}| < 0.9
 * R_{out}, \\ \frac{R_{in}(R_{out} - |\vec{\xi}|)(dE_{a}(t) -
 * dE_{b}(t))}{|\vec{\xi}|(R_{out} - R_{in})}, & 0.9 R_{out} \leq |\vec{\xi}| <
 * R_{out} \end{array}\right.
 * \f}
 *
 * ## Jacobian
 * The jacobian is also found through different equations based on which maps
 * are supplied.
 *
 * ### If Rotation, Expansion and Translation maps are supplied then the
 * jacobian is found through
 *
 * \f{equation}{
 * {J^{i}}_{j} = \left\{\begin{array}{ll}E_{a}(t){R^{i}}_{j}(t), & |\vec{\xi}|
 * \leq R_{in},
 * \\ {R^{i}}_{j}(t)E_{a}(t) + \frac{\alpha
 * {R^{i}}_{l}(t)\vec{\xi}^{l}\vec{\xi}_{j}(E_{A}(t) - E_{B}(t))}{\vec{\xi}^2} +
 * w_{E}{R^{i}}_{j}(t), & R_{in} < |\vec{\xi}| < 0.9 R_{out},
 * \\ {R^{i}}_{j}(t)E_{b}(t) + \frac{\alpha
 * {R^{i}}_{l}(t)\vec{\xi}^{l}\vec{\xi}_{j}(E_{A}(t) - E_{B}(t))}{\vec{\xi}^2} +
 * w_{E}{R^{i}}_{j}(t), & 0.9 R_{out} \leq |\vec{\xi}| < R_{out},
 * \\ E_{a}(t){R^{i}}_{j}(t), & |\vec{\xi}| \geq R_{out} \end{array}\right.
 * \f}
 *
 * where $\alpha = \frac{R_{in}R_{out}}{\vec{\xi}^2(R_{in} - R_{out})}$
 *
 * \note For the translation map, the map returns the identity for all regions
 * except between $R_{in}$ and $R_{out}$
 *
 * ## Inverse Jacobian
 * The inverse jacobian is computed numerically by inverting the jacobian.
 */
template <size_t Dim>
class RotScaleTrans {
 public:
  static constexpr size_t dim = Dim;

  RotScaleTrans() = default;

  explicit RotScaleTrans(std::pair<std::string, std::string> scale_f_of_t_names,
                         std::string rot_f_of_t_name,
                         std::string trans_f_of_t_name, double inner_radius,
                         double outer_radius, bool rigid_translation);
  explicit RotScaleTrans(std::pair<std::string, std::string> scale_f_of_t_names,
                         double inner_radius, double outer_radius);
  explicit RotScaleTrans(std::string rot_f_of_t_name);
  explicit RotScaleTrans(std::string trans_f_of_t_name, double inner_radius,
                         double outer_radius, bool rigid_translation);
  explicit RotScaleTrans(std::pair<std::string, std::string> scale_f_of_t_names,
                         std::string rot_f_of_t_name, double inner_radius,
                         double outer_radius);
  explicit RotScaleTrans(std::string rot_f_of_t_name,
                         std::string trans_f_of_t_name, double inner_radius,
                         double outer_radius, bool rigid_translation);
  explicit RotScaleTrans(std::pair<std::string, std::string> scale_f_of_t_names,
                         std::string trans_f_of_t_name, double inner_radius,
                         double outer_radius, bool rigid_translation);

  RotScaleTrans(const RotScaleTrans<Dim>& RotScaleTrans_Map);

  ~RotScaleTrans() = default;
  RotScaleTrans(RotScaleTrans&&) = default;
  RotScaleTrans& operator=(RotScaleTrans&&) = default;
  RotScaleTrans& operator=(const RotScaleTrans& RotScaleTrans_Map);

  template <typename T>
  std::array<tt::remove_cvref_wrap_t<T>, Dim> operator()(
      const std::array<T, Dim>& source_coords, double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) const;

  /// The inverse function is only callable with doubles because the inverse
  /// might fail if called for a point out of range, and it is unclear
  /// what should happen if the inverse were to succeed for some points in a
  /// DataVector but fail for other points.
  std::optional<std::array<double, Dim>> inverse(
      const std::array<double, Dim>& target_coords, double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) const;

  template <typename T>
  std::array<tt::remove_cvref_wrap_t<T>, Dim> frame_velocity(
      const std::array<T, Dim>& source_coords, double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) const;

  template <typename T>
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame> inv_jacobian(
      const std::array<T, Dim>& source_coords, double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) const;

  template <typename T>
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame> jacobian(
      const std::array<T, Dim>& source_coords, double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  static bool is_identity() { return false; }

  const std::unordered_set<std::string>& function_of_time_names() const {
    return f_of_t_names_;
  }

 private:
  template <size_t LocalDim>
  friend bool operator==(  // NOLINT(readability-redundant-declaration)
      const RotScaleTrans<LocalDim>& lhs, const RotScaleTrans<LocalDim>& rhs);

  double root_helper(std::optional<std::array<double, 2>> roots) const;

  std::optional<std::string> scale_f_of_t_a_{};
  std::optional<std::string> scale_f_of_t_b_{};
  std::optional<std::string> rot_f_of_t_{};
  std::optional<std::string> trans_f_of_t_{};
  std::unordered_set<std::string> f_of_t_names_;
  std::optional<double> inner_radius_;
  std::optional<double> outer_radius_;
  bool rigid_ = true;
};

template <size_t Dim>
bool operator!=(const RotScaleTrans<Dim>& lhs, const RotScaleTrans<Dim>& rhs) {
  return not(lhs == rhs);
}

}  // namespace domain::CoordinateMaps::TimeDependent
