// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/CoordinateMaps/Distribution.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Utilities/TypeTraits/RemoveReferenceWrapper.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace domain::CoordinateMaps {

namespace detail {
// This mapping can be deleted once the 2D and 3D wedges are oriented the same
// (see issue https://github.com/sxs-collaboration/spectre/issues/2988)
template <size_t Dim>
struct WedgeCoordOrientation;
template <>
struct WedgeCoordOrientation<2> {
  static constexpr size_t radial_coord = 0;
  static constexpr size_t polar_coord = 1;
  static constexpr size_t azimuth_coord = 2;  // unused
};
template <>
struct WedgeCoordOrientation<3> {
  static constexpr size_t radial_coord = 2;
  static constexpr size_t polar_coord = 0;
  static constexpr size_t azimuth_coord = 1;
};
}  // namespace detail

// TODO : do we need to assert anything regarding cube_half_length_ and
// focal_offset_? or cube_half_length_ and radius_inner_ or radius_outer_?
// TODO : make sure everything still makes sense for variable opening angles
// even with an offset (do in review with Marcie)
// TODO : add pictures of centered and offset wedge with things labeled
// TODO : run tests in Debug
// TODO : document jacobian and inverse jacobian with offset, include gamma
// TODO: ask Marcie if we need to disambiguate the two different cap Xi and
// cap eta
// TODO : do we need to make a point to explain that when opening angles are
// pi/2, whether with_adapted_equiangular_map_ is true or false changes
// nothing? maybe in the new opening angles section where it doscusses
// theta_D ?
// TODO : make sure cling-tidy and CI pass
/*!
 * \ingroup CoordinateMapsGroup
 *
 * \brief Map from a square or cube to a wedge.
 * \image html Shell.png "A shell can be constructed out of six wedges."
 *
 * \details The mapping that goes from a reference cube (in 3D) or square (in
 *  2D) to a wedge centered on a coordinate axis covering a volume between an
 *  inner surface and outer surface. Each surface can be given a curvature
 *  between flat (a sphericity of 0) or spherical (a sphericity of 1).
 *
 *  In 2D, the first logical coordinate corresponds to the radial coordinate,
 *  and the second logical coordinate corresponds to the angular coordinate. In
 *  3D, the first two logical coordinates correspond to the two angular
 *  coordinates, and the third to the radial coordinate. This difference
 *  originates from separate implementations for the 2D and 3D map that were
 *  merged. The 3D implementation can be changed to use the first logical
 *  coordinate as the radial direction, but this requires propagating the change
 *  through the rest of the domain code (see issue
 *  https://github.com/sxs-collaboration/spectre/issues/2988).
 *
 *  The following documentation is for the **centered** 3D map, as we will defer
 *  the dicussion of `Wedge`s with a `focal_offset_` to a later section. The 2D
 *  map is obtained by setting either of the two angular coordinates to zero
 *  (and using \f$\xi\f$ as the radial coordinate).
 *
 *  The Wedge map is constructed by linearly interpolating between a bulged
 *  face of radius `radius_inner_` to a bulged face of radius `radius_outer_`,
 *  where the radius of each bulged face is defined to be the radius of the
 *  sphere circumscribing the bulge.
 *
 *  We make a choice here as to whether we wish to use the logical coordinates
 *  parameterizing these surface as they are, in which case we have the
 *  equidistant choice of coordinates, or whether to apply a tangent map to them
 *  which leads us to the equiangular choice of coordinates. `Wedge`s have
 *  variable `opening_angles_` which, for centered `Wedge`s, are the angular
 *  sizes of the wedge in the $\xi$ and $\eta$ directions (for the 3D case) in
 * the target frame. By default, `Wedge`s have opening angles of $\pi/2$, so we
 * will discuss that case here and defer both the discussion of generalized
 * opening angles and the interaction between opening angles and non-zero focal
 * offsets for later sections.
 *
 *  For a Wedge with $\xi$ and $\eta$ opening angles of $\pi/2$, the
 *  equiangular coordinates in terms of the logical coordinates are:
 *
 *  \begin{align}
 *    \textrm{equiangular xi} : \Xi(\xi) = \textrm{tan}(\xi\pi/4)
 *    \label{eq:equiangular_xi_pi_over_2}
 *  \end{align}
 *
 *  \begin{align}
 *    \textrm{equiangular eta} :
 *        \mathrm{H}(\eta) =  \textrm{tan}(\eta\pi/4)
 *        \label{eq:equiangular_eta_pi_over_2}
 *  \end{align}
 *
 *  With derivatives:
 *
 *  \begin{align}
 *    \Xi'(\xi) &= \frac{\pi}{4}(1+\Xi^2) \\
 *    \mathrm{H}'(\eta) &= \frac{\pi}{4}(1+\mathrm{H}^2)
 *  \end{align}
 *
 *  The equidistant coordinates are:
 *
 *  \begin{align}
 *    \textrm{equidistant xi}  : \Xi = \xi \\
 *    \textrm{equidistant eta}  : \mathrm{H} = \eta
 *  \end{align}
 *
 *  with derivatives:
 *
 *  \begin{align}
 *    \Xi'(\xi) &= 1 \\
 *    \mathrm{H}'(\eta) &= 1
 *  \end{align}
 *
 *  We also define the variable \f$\rho\f$, given by:
 *
 *  \begin{align}
 *    \textrm{rho} : \rho = \sqrt{1+\Xi^2+\mathrm{H}^2}
 *  \end{align}
 *
 *  ### The Spherical Face Map
 *  The surface map for the spherical face of radius \f$R\f$ lying in the
 *  \f$+z\f$ direction in either choice of coordinates is then given by:
 *
 *  \begin{align}
 *    \vec{\sigma}_{spherical}: \vec{\xi} \rightarrow \vec{x}(\vec{\xi})
 *  \end{align}
 *
 *  Where
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta) =
 *        \begin{bmatrix}
 *          x(\xi,\eta) \\
 *          y(\xi,\eta) \\
 *          z(\xi,\eta) \\
 *        \end{bmatrix}  =
 *            \frac{R}{\rho}
 *                \begin{bmatrix}
 *                  \Xi \\
 *                  \mathrm{H} \\
 *                  1 \\
 *                \end{bmatrix}
 *  \end{align}
 *
 *  ### The Bulged Face Map
 *  The bulged surface is itself constructed by linearly interpolating between
 *  a cubical face and a spherical face. The surface map for the cubical face
 *  of side length \f$2L\f$ lying in the \f$+z\f$ direction is given by:
 *
 *  \begin{align}
 *    \vec{\sigma}_{cubical}: \vec{\xi} \rightarrow \vec{x}(\vec{\xi})
 *  \end{align}
 *
 *  Where
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta) =
 *        \begin{bmatrix}
 *          x(\xi,\eta) \\
 *          y(\xi,\eta) \\
 *          L \\
 *        \end{bmatrix} =
 *            L\begin{bmatrix}
 *               \Xi \\
 *               \mathrm{H} \\
 *               1 \\
 *             \end{bmatrix}
 *  \end{align}
 *
 *  To construct the bulged map we interpolate between this cubical face map
 *  and a spherical face map of radius \f$R\f$, with the interpolation
 *  parameter being \f$s\f$, called the *sphericity* and which ranges from
 *  0 to 1, with 0 corresponding to a flat surface and 1 corresponding to a
 *  spherical surface. The surface map for the bulged face lying in the \f$+z\f$
 *  direction is then given by:
 *
 *  \begin{align}
 *    \vec{\sigma}_{bulged}(\xi,\eta) =
 *        \left\{(1-s)L +
 *        \frac{sR}{\rho}\right\}
 *            \begin{bmatrix}
 *              \Xi \\
 *              \mathrm{H} \\
 *              1 \\
 *            \end{bmatrix}
 *  \end{align}
 *
 *  We constrain $L$ by demanding that the spherical face circumscribe the cube.
 *  With this condition, we have \f$L = R/\sqrt3\f$.
 *  \note This differs from the choice in SpEC where it is demanded that the
 *  surfaces touch at the center, which leads to \f$L = R\f$.
 *
 *  ### The Full Volume Map
 *  The final map for the wedge which lies along the \f$+z\f$ axis is obtained
 *  by interpolating between the two surfaces with the interpolation parameter
 *  being the logical coordinate \f$\zeta\f$. For a wedge whose gridpoints are
 *  **linearly** distributed in the radial direction (`radial_distribution_` is
 *  \ref domain::CoordinateMaps::Distribution
 *  "domain::CoordinateMaps::Distribution::Linear"), this interpolation results
 *  in the following map:
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta,\zeta) =
 *        \frac{1}{2}\left\{
 *          (1-\zeta)\Big[
 *            (1-s_{inner})\frac{R_{inner}}{\sqrt 3} +
 *            s_{inner}\frac{R_{inner}}{\rho}
 *          \Big] +
 *          (1+\zeta)\Big[
 *            (1-s_{outer})\frac{R_{outer}}{\sqrt 3} +
 *            s_{outer}\frac{R_{outer}}{\rho}
 *          \Big]
 *        \right\}
 *            \begin{bmatrix}
 *              \Xi \\
 *              \mathrm{H} \\
 *              1 \\
 *            \end{bmatrix}
 *  \end{align}
 *
 *  We will define the variables \f$F(\zeta)\f$ and \f$S(\zeta)\f$, the frustum
 *  and sphere factors (in the linear case):
 *
 *  \begin{align}
 *    F(\zeta) &= F_0 + F_1\zeta \\
 *    S(\zeta) &= S_0 + S_1\zeta
 *  \end{align}
 *
 *  Where
 *
 *  \begin{align}
 *    F_0 &=
 *        \frac{1}{2} \big\{
 *          (1-s_{outer})R_{outer} + (1-s_{inner})R_{inner}
 *        \big\} \\
 *    F_1 &= \partial_{\zeta}F
 *         = \frac{1}{2} \big\{
 *             (1-s_{outer})R_{outer} - (1-s_{inner})R_{inner}
 *           \big\} \\
 *    S_0 &=
 *        \frac{1}{2} \big\{
 *          s_{outer}R_{outer} + s_{inner}R_{inner}
 *        \big\} \\
 *    S_1 &= \partial_{\zeta}S
 *         = \frac{1}{2} \big\{ s_{outer}R_{outer} - s_{inner}R_{inner}\big\}
 *  \end{align}
 *
 *  The map can then be rewritten as:
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta,\zeta) =
 *        \left\{
 *          \frac{F(\zeta)}{\sqrt 3} + \frac{S(\zeta)}{\rho}
 *        \right\}
 *            \begin{bmatrix}
 *              \Xi \\
 *              \mathrm{H} \\
 *              1 \\
 *            \end{bmatrix}
 *  \end{align}
 *
 *  The inverse map is given by:
 *
 *  \begin{align}
 *    \xi &= \frac{x}{z} \\
 *    \eta &= \frac{y}{z} \\
 *    \zeta &= \frac{z - \left(\frac{S_0}{\rho} + \frac{F_0}{\sqrt{3}}\right)}
 *                  {\left(\frac{S_1}{\rho} + \frac{F_1}{\sqrt{3}}\right)}
 *  \end{align}
 *
 *  We provide some common derivatives:
 *
 *  \begin{align}
 *    \partial_{\xi}z &= \frac{-S(\zeta)\Xi\Xi'}{\rho^3} \\
 *    \partial_{\eta}z &= \frac{-S(\zeta)\mathrm{H}\mathrm{H}'}{\rho^3} \\
 *    \partial_{\zeta}z &= \frac{F'}{\sqrt 3} + \frac{S'(\zeta)}{\rho}
 *  \end{align}
 *
 *  The Jacobian then is:
 *
 *  \begin{align}
 *    J =
 *        \begin{bmatrix}
 *          \Xi'z + \Xi\partial_{\xi}z &
 *              \Xi\partial_{\eta}z &
 *              \Xi\partial_{\zeta}z \\
 *          \mathrm{H}\partial_{\xi}z &
 *              \mathrm{H}'z + \mathrm{H}\partial_{\eta}z &
 *              \mathrm{H}\partial_{\zeta}z \\
 *          \partial_{\xi}z &
 *              \partial_{\eta}z &
 *              \partial_{\zeta}z \\
 *        \end{bmatrix}
 *  \end{align}
 *
 *  A common factor that shows up in the inverse jacobian is:
 *  \begin{align}
 *    T:= \frac{S(\zeta)}{(\partial_{\zeta}z)\rho^3}
 *  \end{align}
 *
 *  The inverse Jacobian then is:
 *  \begin{align}
 *    J^{-1} =
 *        \frac{1}{z}\begin{bmatrix}
 *          \Xi'^{-1} & 0 & -\Xi\Xi'^{-1} \\
 *          0 & \mathrm{H}'^{-1} & -\mathrm{H}\mathrm{H}'^{-1} \\
 *          T\Xi & T\mathrm{H} & T + F(\partial_{\zeta}z)^{-1}/\sqrt 3 \\
 *        \end{bmatrix}
 *  \end{align}
 *
 *  ### Changing the radial distribution of the gridpoints
 *  By default, Wedge linearly distributes its gridpoints in the radial
 *  direction. An exponential distribution of gridpoints can be obtained by
 *  linearly interpolating in the logarithm of the radius in order to obtain
 *  a relatively higher resolution at smaller radii. Since this is a radial
 *  rescaling of Wedge, this option is only supported for fully spherical
 *  wedges with `sphericity_inner_` = `sphericity_outer_` = 1.
 *
 *  The linear interpolation done for a logarithmic radial distribution
 *  (`radial_distribution_` is \ref domain::CoordinateMaps::Distribution
 *  "domain::CoordinateMaps::Distribution::Logarithmic") is:
 *
 *  \begin{align}
 *    \ln r = \frac{1-\zeta}{2}\ln R_{inner} + \frac{1+\zeta}{2}\ln R_{outer}
 *  \end{align}
 *
 *  The map then is:
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta,\zeta) =
 *        \frac{\sqrt{R_{inner}^{1-\zeta}R_{outer}^{1+\zeta}}}{\rho}
 *            \begin{bmatrix}
 *              \Xi \\
 *              \mathrm{H} \\
 *              1 \\
 *            \end{bmatrix}
 *  \end{align}
 *
 *  We can rewrite this map to take on the same form as the map for the linear
 *  radial distribution, where we set
 *
 *  \begin{align}
 *    F(\zeta) &= 0 \\
 *    S(\zeta) &= \sqrt{R_{inner}^{1-\zeta}R_{outer}^{1+\zeta}} \\
 *  \end{align}
 *
 *  Which gives us
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta,\zeta) =
 *        \frac{S(\zeta)}{\rho}
 *            \begin{bmatrix}
 *              \Xi \\
 *              \mathrm{H} \\
 *              1 \\
 *            \end{bmatrix}
 *  \end{align}
 *
 *  TODO: Reference equation
 *  The jacobian then also takes the same form in terms of this \f$S(\zeta)\f$
 *  and \f$S'(\zeta)\f$.
 *
 *  Alternatively, an inverse radial distribution (`radial_distribution_` is
 *  \ref domain::CoordinateMaps::Distribution
 *  "domain::CoordinateMaps::Distribution::Inverse") can be chosen where the
 *  linear interpolation is:
 *
 *  \begin{align}
 *    \frac{1}{r} =
 *        \frac{R_\mathrm{inner} + R_\mathrm{outer}}
 *             {2 R_\mathrm{inner}R_\mathrm{outer}} +
 *        \frac{R_\mathrm{inner} - R_\mathrm{outer}}
 *             {2R_\mathrm{inner} R_\mathrm{outer}} \zeta
 *  \end{align}
 *
 *  Which can be rewritten as:
 *
 *  \begin{align}
 *    \frac{1}{r} = \frac{1-\zeta}{2R_{inner}} + \frac{1+\zeta}{2R_{outer}}
 *  \end{align}
 *
 *  The map likewise takes the form:
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta,\zeta) =
 *        \frac{S(\zeta)}{\rho}
 *            \begin{bmatrix}
 *              \Xi \\
 *              \mathrm{H} \\
 *              1 \\
 *            \end{bmatrix}
 *  \end{align}
 *
 *  Where
 *
 *  \begin{align}
 *    S(\zeta) =
 *        \frac{2R_{inner}R_{outer}}
 *             {(1 + \zeta)R_{inner} + (1 - \zeta)R_{outer}}
 *  \end{align}
 *
 *  TODO: Reference equation
 *  And the jacobian again takes the same form in terms of this \f$S(\zeta)\f$
 *  and \f$S'(\zeta)\f$.
 *
 *  ### Changing the opening angles
 *  Consider the following map on $\xi \in [-1,1]$, which maps this interval
 *  onto a parameterized curve that extends one fourth of a circle.
 *
 *  \begin{align}
 *    \vec{\Gamma}(\xi) =
 *        \frac{R}{\sqrt{1+\xi^2}}
 *            \begin{bmatrix}
 *              1 \\
 *              \xi \\
 *              \end{bmatrix}.
 *    \label{eq:quarter_circle}
 *  \end{align}
 *
 *  It is convenient to compute the polar coordinate $\theta$ of the mapped
 *  point as a function of $\xi$:
 *
 *  \begin{align}
 *    \theta(\xi) = \tan^{-1}\left(\frac{\Gamma_y(\xi)}{\Gamma_x(\xi)}\right).
 *    \label{eq:polar_coord}
 *  \end{align}
 *
 * The *opening angle* of the map is defined to be:
 *
 *  \begin{align}
 *    \Delta \theta = \theta(1) - \theta(-1),
 *    \label{eq:define_opening_angle}
 *  \end{align}
 *
 *  We can see that with $\xi=\pm 1$, we have $\Gamma_x = R/\sqrt{2}$ and
 *  $\Gamma_y=\pm R/\sqrt{2}$, giving us
 *  $\theta(1) = \pi/4$ and $\theta(-1) = -\pi/4$. This wedge has an opening
 *  angle $\pi/2$ radians, as expected. On the other hand, the following map has
 *  an opening angle of $\theta_O$:
 *
 *  \begin{align}
 *    \vec{\Gamma}(\xi) =
 *        \frac{R}{\sqrt{1+\tan^2{(\theta_O/2)}\xi^2}}
 *            \begin{bmatrix}
 *            1 \\
 *            \tan{(\theta_O/2)}\xi \\
 *            \end{bmatrix}.
 *  \end{align}
 *
 *  Let us also consider the generalized map
 *
 *  \begin{align}
 *    \vec{\Gamma}(\xi) =
 *        \frac{R}{\sqrt{1+\Xi^2}}
 *            \begin{bmatrix}
 *              1 \\
 *              \Xi \\
 *            \end{bmatrix},
 *  \end{align}
 *
 *  where $\Xi(\xi)$ is a function of $\xi$. $\theta(\xi)$ can then be written
 *  as
 *
 *  \begin{align}
 *    \theta(\xi) = \tan^{-1}(\Xi).
 *  \end{align}
 *
 *  A curve $\vec{\Gamma}(\xi)$ is parameterized equiangularly if
 *
 *  \begin{align}
 *    \frac{\mathrm{d}\theta}{\mathrm{d}\xi} = C.
 *  \end{align}
 *
 *  For the equiangular map with an opening angle of $pi/2$,
 *  $\Xi = \tan(\pi\xi/4)$, $\theta = \pi\xi/4$, and $\Delta\theta = \pi/2$. As
 *  for the map
 *
 *  \begin{align}
 *    \Xi(\xi) =
 *        \tan{(\theta_O/2)}\frac{\tan{(\theta_D \xi/2)}}{\tan{(\theta_D/2)}},
 *  \end{align}
 *
 *  this choice of $\Xi(\xi)$ results in a $\vec{\Gamma}(\xi)$ with opening
 *  angle $\theta_O$, which is equiangularly distributed if
 *  $\theta_O = \theta_D$. In the Wedge map, the argument
 *  `with_adapted_equiangular_map` controls whether to set
 *  $\theta_O = \theta_D$ (the `true` case) or to set $\theta_D = \pi/2$
 *  (the `false` case). When working with a 3D Wedge, the opening angles for the
 *  Wedge can be separately controlled for both the $\xi$ and $\eta$ directions,
 *  but `with_adapted_equiangular_map` will apply to both directions.
 *  Additionally in the 3D case, it is not possible to set
 *  `with_equiangular_map_` to `true` for all of the six wedges of a sphere
 *  unless every opening angle is $\pi/2$. In the
 *  \ref ::domain::creators::BinaryCompactObject "BinaryCompactObject" domain,
 *  the outer $+y$, $-y$, $+z$, and $-z$ `Wedge`s are allowed to have a
 *  user-specified opening angle in the $\xi$-direction, with a corresponding
 *  $\theta_D$ equal to this opening angle, while in the $\eta$-direction the
 *  opening angle is set to $\pi/2$. The two end cap `Wedge`s in the $+x$ and
 *  $-x$ directions have angular dimensions and gridpoint distributions
 *  determined by the other four `Wedge`s, as the six `Wedge`s must conforming
 *  have gridpoint distributions at the $\xi = \pm1$, $\eta = \pm 1$ boundaries.
 *
 *  ### Wedge with a Focal Offset
 *  In the case of the rectangular
 *  \ref ::domain::creators::BinaryCompactObject "BinaryCompactObject" domain,
 *  it becomes desirable to offset the center of the spherical excision surface
 *  relative to the center of the cubical surface surrounding it. To enable the
 *  offsetting of the central excision, the Wedge map must be generalized
 *  according to the *focal lifting* method, which we will now discuss.
 *
 *  We consider the problem of creating parameterized volumes from parameterized
 *  surfaces. Consider a parameterized surface $\vec{\rho}(\xi,\eta)$, also
 *  referred to as the *parent surface*. We define *focal lifting* as the
 *  projection of this parent surface into a three-dimensional parameterized
 *  volume $\vec{x}(\xi,\eta, \zeta)$ with respect to some *focus* $\vec{x}_0$
 *  and *lifting scale factor* $\Lambda(\xi,\eta,\zeta)$. The resulting volume
 *  is then said to be a *focally lifted* volume. These volume maps can be cast
 *  into the following form:
 *
 *  \begin{align}
 *    \vec{x} - \vec{x}_0 = \Lambda(\vec{\rho}-\vec{x}_0),
 *  \end{align}
 *
 *  which makes apparent how the mapped point $\vec{x}(\xi,\eta,\zeta)$ is
 *  obtained. The parametric equations for the generalized 3D Wedge maps can all
 *  be written in the above form, which we will refer to as
 *  *focally lifted form*. In the case of the 3D Wedge map with no focal offset,
 *  we have:
 *
 *  \begin{align}
 *    \vec{x}_0 &= 0 \\
 *    \Lambda &= \left\{\frac{F(\zeta)}{\sqrt{3}} +
 *                      \frac{S(\zeta)}{\rho} \right\} \\
 *    \vec{\rho} &= \begin{bmatrix} \Xi, \mathrm{H}, 1 \end{bmatrix}^T
 *  \end{align}
 *
 *  The above map can be thought of as constructing a wedge from a biunit cube
 *  centered at the origin. Points on the parent surface are scaled by a factor
 *  of $\Lambda(\xi,\eta,\zeta)$ to obtain the corresponding point in the
 *  volume. When generalizing the map to have a focus shifted from the origin
 *  (obtained by setting `focal_offset_` to be non-zero), we scale the original
 *  parent surface $\vec{\rho} = [\Xi, \mathrm{H},1]^T$ by a factor $L$, and let
 *  the focus $\vec{x_0}$ shift away from the origin. The generalized wedge map
 *  is then given by:
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta,\zeta) =
 *        \left\{\frac{F(\zeta)}{L\sqrt 3} +
 *        \frac{S(\zeta)}{L\rho}\right\}
 *            \begin{bmatrix}
 *              L\Xi - x_0 \\
 *              L\mathrm{H} - y_0 \\
 *              L-z_0 \\
 *            \end{bmatrix}
 *  \end{align}
 *
 *  where $\rho$ is now
 *  $\sqrt{(\Xi - x_0/L)^2 + (\mathrm{H} - y_0/L)^2 + (1 - z_0/L)^2}$.
 *
 *  This map is often written as:
 *
 *  \begin{align}
 *    \vec{x}(\xi,\eta,\zeta) =
 *        \left\{\frac{F(\zeta)}{\sqrt{3}} +
 *        \frac{S(\zeta)}{\rho}\right\}(\vec{\sigma}_0 - \vec{x}_0/L),
 *    \label{eq:focally_lifted_map_with_generalized_z_coef}
 *  \end{align}
 *
 *  where $\vec{\sigma}_0 = [\Xi, \mathrm{H},1]^T$, as the parent surface
 *  $\vec{\rho}$ is now $L\vec{\sigma}_0$. We give the quantity in braces the
 *  name $z_{\Lambda} = L\Lambda$, *generalized z*. The map can be inverted by
 *  first solving for \f$z_{\Lambda}\f$ in terms of the target coordinates. We
 *  make use of the fact that the parent surface $\vec{\rho}$ has a constant
 *  normal vector $\hat{n} = \hat{z}$.
 *
 *  \begin{align}
 *    z_{\Lambda} = \frac{(\vec{x} - \vec{x}_0)\cdot\hat{n}}
 *                       {(\vec{\sigma}_0-\vec{x}_0/L)\cdot\hat{n}}.
 *  \end{align}
 *
 *  Moving all the known quantities to the left hand side results in the
 *  following:
 *
 *  \begin{align}
 *    \frac{\vec{x} - \vec{x}_0}{z_{\Lambda}} + \frac{\vec{x}_0}{L}
 *         = \vec{\sigma}_0(\xi,\eta)
 *         = \begin{bmatrix}
 *             \Xi \\
 *             \mathrm{H} \\
 *             1 \\
 *           \end{bmatrix},
 *  \end{align}
 *
 *  Note that $|\vec{\sigma}_0 - \vec{x}_0/L| = \sqrt{(\Xi - x_0/L)^2 +
 *  (\mathrm{H} - y_0/L)^2 + (1 - z_0/L)^2} = \rho$, indicating that an
 *  expression for $\rho$ in terms of the target coordinates can be computed via
 *  taking the magnitude of both sides of
 *  Eq. ($\ref{eq:focally_lifted_map_with_generalized_z_coef}$):
 *
 *  \begin{align}
 *    |\vec{x} - \vec{x}_0| = z_{\Lambda}|\vec{\sigma}_0 - \vec{x}_0/L|
 *                          = z_{\Lambda}\rho.
 *  \end{align}
 *
 *  The quantity $\rho$ is then given by:
 *
 *  \begin{align}
 *    \rho = \frac{|\vec{x} - \vec{x}_0|}{z_{\Lambda}}.
 *  \end{align}
 *
 *  With $\rho$ computed, $\zeta$ can be computed from
 *
 *  \begin{align}
 *    z_{\Lambda} = \left\{\frac{F(\zeta)}{\sqrt{3}} +
 *                  \frac{S(\zeta)}{\rho} \right\}
 *                = \left\{\frac{F_0}{\sqrt{3}} + \frac{S_0}{\rho} +
 *                  \frac{F_1\zeta}{\sqrt{3}} + \frac{S_1\zeta}{\rho}\right\},
 *  \end{align}
 *
 *  which gives
 *
 *  \begin{align}
 *     \zeta = \frac{z_{\Lambda} -
 *                   \left(\frac{S_0}{\rho} + \frac{F_0}{\sqrt{3}}\right)}
 *                  {\left(\frac{S_1}{\rho} + \frac{F_1}{\sqrt{3}}\right)}.
 *  \end{align}
 *
 *  TODO : jacobian stuff here
 *
 *  ### Interaction between opening angles and focal offsets
 *  When a Wedge is created with a non-zero focal offset, the resulting shape
 *  can take on a variety of possible angular sizes, depending on where the
 *  focus is placed relative to the default centered location. The reader might
 *  note that the angular size of a Wedge can also be controlled by passing an
 *  argument to the `opening_angles` parameter in the Wedge constructor. While
 *  both of these methods allow the angular size of a Wedge to be changed, the
 *  user is prevented from employing both of them at the same time. In
 *  particular, when the the offset is set to some non-zero value, the
 *  `opening_angles_` member variable is set to $\pi/2$. Note that the
 *  `opening_angles_` member being set to $\pi/2$ does not imply the
 *  resulting Wedge will have an angular size of $\pi/2$. On the contrary, the
 *  Wedge will have the angular size that is determined by the application of
 *  the focal lifting method on the parent surface, which is the upper $+z$ face
 *  of a cube that is centered at the origin.
 *
 *  Because `opening_angles_` is set to $\pi/2$ when there is a non-zero focal
 *  offset, when there is a non-zero focal offset and `with_equiangular_map_` is
 *  `true`, $\Xi$ is given by Eq. ($\ref{eq:equiangular_xi_pi_over_2}$) and
 *  $\mathrm{H}$ by Eq. ($\ref{eq:equiangular_eta_pi_over_2}$), just as it is
 *  for the case of a centered Wedge with `opening_angles_` of $\pi/2$.
 */
template <size_t Dim>
class Wedge {
 public:
  static constexpr size_t dim = Dim;
  enum class WedgeHalves {
    /// Use the entire wedge
    Both,
    /// Use only the upper logical half
    UpperOnly,
    /// Use only the lower logical half
    LowerOnly
  };

  /*!
   * Constructs a 3D wedge.
   * \param radius_inner Distance from the origin to one of the
   * corners which lie on the inner surface.
   * \param radius_outer Distance from the origin to one of the
   * corners which lie on the outer surface.
   * \param orientation_of_wedge The orientation of the desired wedge relative
   * to the orientation of the default wedge which is a wedge that has its
   * curved surfaces pierced by the upper-z axis. The logical xi and eta
   * coordinates point in the cartesian x and y directions, respectively.
   * \param sphericity_inner Value between 0 and 1 which determines
   * whether the inner surface is flat (value of 0), spherical (value of 1) or
   * somewhere in between
   * \param sphericity_outer Value between 0 and 1 which determines
   * whether the outer surface is flat (value of 0), spherical (value of 1) or
   * somewhere in between
   * \param cube_half_length Half the length of the parent surface (see Wedge
   * documentation for more details)
   * \param focal_offset The target frame coordinates of the focus from which
   * the Wedge is focally lifted
   * \param with_equiangular_map Determines whether to apply a tangent function
   * mapping to the logical coordinates (for `true`) or not (for `false`).
   * \param halves_to_use Determines whether to construct a full wedge or only
   * half a wedge. If constructing only half a wedge, the resulting shape has a
   * face normal to the x direction (assuming default OrientationMap). If
   * constructing half a wedge, an intermediate affine map is applied to the
   * logical xi coordinate such that the interval [-1,1] is mapped to the
   * corresponding logical half of the wedge. For example, if `UpperOnly` is
   * specified, [-1,1] is mapped to [0,1], and if `LowerOnly` is specified,
   * [-1,1] is mapped to [-1,0]. The case of `Both` means a full wedge, with no
   * intermediate map applied. In all cases, the logical points returned by the
   * inverse map will lie in the range [-1,1] in each dimension. Half wedges are
   * currently only useful in constructing domains for binary systems.
   * \param radial_distribution Determines how to distribute gridpoints along
   * the radial direction. For wedges that are not exactly spherical, only
   * `Distribution::Linear` is currently supported.
   * \param opening_angles Determines the angular size of the wedge. The
   * default value is pi/2, which corresponds to a wedge size of pi/2. For this
   * setting, four Wedges can be put together to cover 2pi in angle along a
   * great circle. This option is meant to be used with the equiangular map
   * option turned on.
   * \param with_adapted_equiangular_map Determines whether to adapt the
   * point distribution in the wedge to match its physical angular size. When
   * `true`, angular distances are proportional to logical distances. Note
   * that it is not possible to use adapted maps in every Wedge of a Sphere
   * unless each Wedge has the same size along both angular directions.
   */
  Wedge(double radius_inner, double radius_outer, double sphericity_inner,
        double sphericity_outer, double cube_half_length,
        std::array<double, Dim> focal_offset,
        OrientationMap<Dim> orientation_of_wedge, bool with_equiangular_map,
        WedgeHalves halves_to_use = WedgeHalves::Both,
        Distribution radial_distribution = Distribution::Linear,
        const std::array<double, Dim - 1>& opening_angles =
            make_array<Dim - 1>(M_PI_2),
        bool with_adapted_equiangular_map = true);

  Wedge() = default;
  ~Wedge() = default;
  Wedge(Wedge&&) = default;
  Wedge(const Wedge&) = default;
  Wedge& operator=(const Wedge&) = default;
  Wedge& operator=(Wedge&&) = default;

  template <typename T>
  std::array<tt::remove_cvref_wrap_t<T>, Dim> operator()(
      const std::array<T, Dim>& source_coords) const;

  /// For a \f$+z\f$-oriented `Wedge`, returns invalid if \f$z<=0\f$
  /// or if \f$(x,y,z)\f$ is on or outside the cone defined
  /// by \f$(x^2/z^2 + y^2/z^2+1)^{1/2} = -S/F\f$, where
  /// \f$S = \frac{1}{2}(s_1 r_1 - s_0 r_0)\f$ and
  /// \f$F = \frac{1}{2\sqrt{3}}((1-s_1) r_1 - (1-s_0) r_0)\f$.
  /// Here \f$s_0,s_1\f$ and \f$r_0,r_1\f$ are the specified sphericities
  /// and radii of the inner and outer \f$z\f$ surfaces.  The map is singular on
  /// the cone and on the xy plane.
  /// The inverse function is only callable with doubles because the inverse
  /// might fail if called for a point out of range, and it is unclear
  /// what should happen if the inverse were to succeed for some points in a
  /// DataVector but fail for other points.
  std::optional<std::array<double, Dim>> inverse(
      const std::array<double, Dim>& target_coords) const;

  template <typename T>
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame> jacobian(
      const std::array<T, Dim>& source_coords) const;

  template <typename T>
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame> inv_jacobian(
      const std::array<T, Dim>& source_coords) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  static constexpr bool is_identity() { return false; }

 private:
  // maps between 2D and 3D choices for coordinate axis orientations
  static constexpr size_t radial_coord =
      detail::WedgeCoordOrientation<Dim>::radial_coord;
  static constexpr size_t polar_coord =
      detail::WedgeCoordOrientation<Dim>::polar_coord;
  static constexpr size_t azimuth_coord =
      detail::WedgeCoordOrientation<Dim>::azimuth_coord;

  /*!
   * \brief Factors out calculation of $S(\zeta)$ needed for the map and the
   * jacobian.
   *
   * \details The value of $S(\zeta)$ is computed differently for different
   * radial distributions.
   *
   * For a **linear** radial distribution:
   *
   * \f{align*}{
   *   S(\zeta) = S_0 + S_1\zeta
   * \f}
   *
   * where $S_0$ and $S_1$ are defined as
   *
   * \f{align*}{
   *   S_0 &=
   *       \frac{1}{2} \big\{
   *         s_{outer}R_{outer} + s_{inner}R_{inner}
   *       \big\} \\
   *   S_1 &= \partial_{\zeta}S
   *        = \frac{1}{2} \big\{ s_{outer}R_{outer} - s_{inner}R_{inner}\big\}
   * \f}
   *
   * and are stored in `sphere_zero_` and `sphere_rate_`, respectively.
   *
   * For a **logarithmic** radial distribution:
   *
   * \f{align*}{
   *   S(\zeta) = \exp{(S_0 + S_1\zeta)}
   * \f}
   *
   * where $S_0$ and $S_1$ are defined as
   *
   * \f{align*}{
   *   S_0 &= \frac{1}{2} \ln(R_{outer}R_{inner}) \\
   *   S_1 &= \frac{1}{2} \ln(R_{outer}/R_{inner})
   * \f}
   *
   * With these definitions of $S_0$ and $S_1$, we can rewrite the expression
   * for $S(\zeta)$ as:
   *
   * \f{align*}{
   *   S(\zeta) &= \sqrt{R_{inner}^{1-\zeta}R_{outer}^{1+\zeta}}
   * \f}
   *
   * As with the linear distribution, $S_0$ and $S_1$ are stored in
   * `sphere_zero_` and `sphere_rate_`, respectively.
   *
   * For an **inverse** radial distribution:
   *
   * \f{align*}{
   *   S(\zeta) =
   *       \frac{2R_{inner}R_{outer}}
   *            {(1 + \zeta)R_{inner} + (1 - \zeta)R_{outer}}
   * \f}
   *
   * In this case, `sphere_zero_` and `sphere_rate_` will simply be `NaN`.
   *
   * See Wedge for more details on these quantities.
   *
   * \param zeta the radial source coordinate
   */
  template <typename T>
  tt::remove_cvref_wrap_t<T> get_s_factor(const T& zeta) const;
  /*!
   * \brief Factors out calculation of $S'(\zeta)$ needed for the jacobian.
   *
   * \details The value of $S'(\zeta)$ is computed differently for different
   * radial distributions.
   *
   * For a **linear** radial distribution:
   *
   * \f{align*}{
   *   S'(\zeta) =
   *       \frac{1}{2} \big\{ s_{outer}R_{outer} - s_{inner}R_{inner}\big\}
   * \f}
   *
   * For a **logarithmic** radial distribution:
   *
   * \f{align*}{
   *   S'(\zeta) = \frac{1}{2} S(\zeta)\ln(R_{outer}/R_{inner})
   * \f}
   *
   * where $S(\zeta)$ is defined in `get_s_factor()`.
   *
   * For an **inverse** radial distribution:
   *
   * \f{align*}{
   *   S'(\zeta) =
   *       \frac{2(R_{inner} R_{outer}^2 - R_{inner}^2 R_{outer})}
   *            {(R_{inner} + R_{outer} + \zeta(R_{inner} - R_{outer}))^2}
   * \f}
   *
   * See Wedge and `get_s_factor()` for more details on these quantities.
   *
   * \param zeta the radial source coordinate
   * \param s_factor $S(\zeta)$ (see `get_s_factor()`)
   */
  template <typename T>
  tt::remove_cvref_wrap_t<T> get_s_factor_deriv(const T& zeta,
                                                const T& s_factor) const;

  /*!
   * \brief Factors out calculation of $z_{\Lambda}$ needed for the map and the
   * jacobian.
   *
   * \details The value of $z_{\Lambda}$  is computed differently for different
   * radial distributions.
   *
   * For a **linear** radial distribution:
   *
   * \f{align*}{
   *   z_{\Lambda} = \frac{F(\zeta)}{\sqrt 3} + \frac{S(\zeta)}{\rho}
   * \f}
   *
   * For a **logarithmic** or **inverse** radial distribution:
   *
   * \f{align*}{
   *   z_{\Lambda} = \frac{S(\zeta)}{\rho}
   * \f}
   *
   * See Wedge and `get_s_factor()` for more details on these quantities.
   *
   * \param zeta the radial source coordinate
   * \param one_over_rho one over $\rho$ where
   * $\rho = |\vec{\sigma}_0 - \vec{x}_0/L| = \sqrt{(\Xi - x_0/L)^2 +
   * (\mathrm{H} - y_0/L)^2 + (1 - z_0/L)^2}$ (see Wedge)
   * \param s_factor $S(\zeta)$ (see `get_s_factor()`)
   */
  template <typename T>
  tt::remove_cvref_wrap_t<T> get_generalized_z(const T& zeta,
                                               const T& one_over_rho,
                                               const T& s_factor) const;
  template <typename T>
  tt::remove_cvref_wrap_t<T> get_generalized_z(const T& zeta,
                                               const T& one_over_rho) const;
  /*!
   * \brief Factors out calculation of $\partial_i z_{\Lambda}$ needed for the
   * jacobian
   *
   * \details For **all** radial distributions:
   *
   * \f{align*}{
   *   \partial_{\xi} z_{\Lambda} &=
   *       \frac{-S(\zeta)\Xi'\gamma_x}{\rho^3} \\
   *   \partial_{\eta} z_{\Lambda} &=
   *       \frac{-S(\zeta)\mathrm{H}'\gamma_y}{\rho^3} \\
   *   \partial_{\zeta} z_{\Lambda} &=
   *       \frac{F'(\zeta)}{\sqrt 3} + \frac{S'(\zeta)}{\rho}
   * \f}
   *
   * However, $\partial_{\zeta} z_{\Lambda}$ reduces to
   *
   * \f{align*}{
   *   \partial_{\zeta} z_{\Lambda} &= \frac{S'(\zeta)}{\rho}
   * \f}
   *
   * for **logarithmic** and **inverse** radial distributions because
   * $F(\zeta) = 0$.
   *
   * See Wedge and `get_s_factor()` for more details on these quantities.
   *
   * \param zeta the radial source coordinate
   * \param one_over_rho one over $\rho$ where
   * $\rho = |\vec{\sigma}_0 - \vec{x}_0/L| = \sqrt{(\Xi - x_0/L)^2 +
   * (\mathrm{H} - y_0/L)^2 + (1 - z_0/L)^2}$ (see Wedge)
   * \param s_factor $S(\zeta)$ (see `get_s_factor()`)
   * \param cap_deriv $\Xi'$ and $\mathrm{H}'$ (see Wedge)
   * \param gamma $\gamma_i$ (see Wedge)
   */
  template <typename T>
  std::array<tt::remove_cvref_wrap_t<T>, Dim> get_d_generalized_z(
      const T& zeta, const T& one_over_rho, const T& s_factor,
      const std::array<tt::remove_cvref_wrap_t<T>, Dim - 1>& cap_deriv,
      const std::array<tt::remove_cvref_wrap_t<T>, Dim>& gamma) const;

  template <size_t LocalDim>
  // NOLINTNEXTLINE(readability-redundant-declaration)
  friend bool operator==(const Wedge<LocalDim>& lhs,
                         const Wedge<LocalDim>& rhs);
  // TODO: Marcie Revisit Constructor Documentation
  /// Distance from the origin to one of the corners which lie on the inner
  /// surface.
  double radius_inner_{std::numeric_limits<double>::signaling_NaN()};
  /// Distance from the origin to one of the corners which lie on the outer
  /// surface.
  double radius_outer_{std::numeric_limits<double>::signaling_NaN()};
  /// Value between 0 and 1 which determines whether the inner surface is flat
  /// (value of 0), spherical (value of 1) or somewhere in between
  double sphericity_inner_{std::numeric_limits<double>::signaling_NaN()};
  /// Value between 0 and 1 which determines whether the outer surface is flat
  /// (value of 0), spherical (value of 1) or somewhere in between
  double sphericity_outer_{std::numeric_limits<double>::signaling_NaN()};
  /// Half the length of the parent surface (see Wedge documentation for more
  /// details)
  double cube_half_length_{std::numeric_limits<double>::signaling_NaN()};
  /// The target frame coordinates of the focus from which the Wedge is focally
  /// lifted
  std::array<double, Dim> focal_offset_{
      make_array<Dim>(std::numeric_limits<double>::signaling_NaN())};
  /// The orientation of the desired wedge relative to the orientation of the
  /// default wedge which is a wedge that has its curved surfaces pierced by the
  /// upper-z axis. The logical xi and eta coordinates point in the cartesian x
  /// and y directions, respectively.
  OrientationMap<Dim> orientation_of_wedge_{};
  /// Determines whether to apply a tangent function mapping to the logical
  /// coordinates (for `true`) or not (for `false`).
  bool with_equiangular_map_ = false;
  /// Determines whether to construct a full wedge or only half a wedge (see
  /// Wedge documentation for more details)
  WedgeHalves halves_to_use_ = WedgeHalves::Both;
  /// Determines how to distribute gridpoints along the radial direction. For
  /// wedges that are not exactly spherical, only `Distribution::Linear` is
  /// currently supported.
  Distribution radial_distribution_ = Distribution::Linear;
  /// $F_0 / \sqrt{3}$ (see Wedge documentation)
  double scaled_frustum_zero_{std::numeric_limits<double>::signaling_NaN()};
  /// $S_0$ (see Wedge documentation)
  double sphere_zero_{std::numeric_limits<double>::signaling_NaN()};
  /// $F_1 / \sqrt{3}$ (see Wedge documentation)
  double scaled_frustum_rate_{std::numeric_limits<double>::signaling_NaN()};
  /// $S_1$ (see Wedge documentation)
  double sphere_rate_{std::numeric_limits<double>::signaling_NaN()};
  /// $\theta_O$ (see Wedge documentation)
  std::array<double, Dim - 1> opening_angles_{
      make_array<Dim - 1>(std::numeric_limits<double>::signaling_NaN())};
  /// $\theta_D$ (see Wedge documentation)
  std::array<double, Dim - 1> opening_angles_distribution_{
      make_array<Dim - 1>(std::numeric_limits<double>::signaling_NaN())};
};

template <size_t Dim>
bool operator!=(const Wedge<Dim>& lhs, const Wedge<Dim>& rhs);
}  // namespace domain::CoordinateMaps
