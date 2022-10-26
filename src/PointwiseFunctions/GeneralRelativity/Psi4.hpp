// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/CoordinateMaps/Tags.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"

/// \cond
namespace gsl {
template <typename>
struct not_null;
}  // namespace gsl
/// \endcond

namespace gr {

/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Computes Newman Penrose quantity $\Psi_4$ using the characteristic
 * field U$^{8+}$ and $\bar{m}^i.
 *
 * \details Computes $\Psi_4$ as: \f$\Psi_4\f$ =
 * U\f$^{8+}_{ij}\bar{m}^i\bar{m}^j\f$ with the characteristic field using the
 * Weyl Electric tensor E_{ij}
 */
template <typename ComplexDataType, size_t SpatialDim, typename Frame,
          typename RealDataType>
void psi_4(
    const gsl::not_null<Scalar<ComplexDataType>*> psi_4_result,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& extrinsic_curvature,
    const tnsr::ijj<RealDataType, SpatialDim, Frame>&
        cov_deriv_extrinsic_curvature,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<RealDataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::I<RealDataType, SpatialDim, Frame>& inertial_coords);

template <typename ComplexDataType, size_t SpatialDim, typename Frame,
          typename RealDataType>
Scalar<ComplexDataType> psi_4(
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& extrinsic_curvature,
    const tnsr::ijj<RealDataType, SpatialDim, Frame>&
        cov_deriv_extrinsic_curvature,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<RealDataType, SpatialDim, Frame>& inverse_spatial_metric,
    // const Scalar<DataType>& sqrt_det_spatial_metric,
    const tnsr::I<RealDataType, SpatialDim, Frame>& inertial_coords);
/// @}
namespace Tags {
/// Computes the Newman Penrose quantity Psi4 using the
/// characteristic field U_{ij}^8+ and m_bar = (x^i - iy^i) / sqrt(2)
///
/// Can be retrieved using `gr::Tags::Psi4`
template <typename ComplexDataType, size_t SpatialDim, typename Frame,
          typename RealDataType>
struct Psi4Compute : Psi4<ComplexDataType>, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::SpatialRicci<SpatialDim, Frame, RealDataType>,
      gr::Tags::ExtrinsicCurvature<SpatialDim, Frame, RealDataType>,
      ::Tags::deriv<gr::Tags::ExtrinsicCurvature<3, Frame, RealDataType>,
                    tmpl::size_t<3>, Frame>,
      gr::Tags::SpatialMetric<SpatialDim, Frame, RealDataType>,
      gr::Tags::InverseSpatialMetric<SpatialDim, Frame, RealDataType>,
      domain::Tags::Coordinates<SpatialDim, Frame>>;

  using return_type = Scalar<ComplexDataType>;
  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<ComplexDataType>*>,
                           const tnsr::ii<RealDataType, SpatialDim, Frame>&,
                           const tnsr::ii<RealDataType, SpatialDim, Frame>&,
                           const tnsr::ijj<RealDataType, SpatialDim, Frame>&,
                           const tnsr::ii<RealDataType, SpatialDim, Frame>&,
                           const tnsr::II<RealDataType, SpatialDim, Frame>&,
                           const tnsr::I<RealDataType, SpatialDim, Frame>&)>(
          &psi_4<ComplexDataType, SpatialDim, Frame, RealDataType>);
  using base = Psi4<ComplexDataType>;
};
}  // namespace Tags
}  // namespace gr
