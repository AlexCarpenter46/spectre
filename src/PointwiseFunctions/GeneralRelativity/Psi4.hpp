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
 *
 *
 */
template <size_t SpatialDim, typename Frame, typename DataType>
void psi_4(
    const gsl::not_null<Scalar<DataType>*> psi_4_result,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<DataType, SpatialDim, Frame>& extrinsic_curvature,
    const tnsr::ijj<DataType, SpatialDim, Frame>& cov_deriv_extrinsic_curvature,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    // const Scalar<DataType>& sqrt_det_spatial_metric,
    const tnsr::I<DataType, SpatialDim, Frame>& inertial_coords);
// const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
// const tnsr::ii<DataType, SpatialDim, Frame>& weyl_magnetic) {

template <size_t SpatialDim, typename Frame, typename DataType>
Scalar<DataType> psi_4(
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<DataType, SpatialDim, Frame>& extrinsic_curvature,
    const tnsr::ijj<DataType, SpatialDim, Frame>& cov_deriv_extrinsic_curvature,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    // const Scalar<DataType>& sqrt_det_spatial_metric,
    const tnsr::I<DataType, SpatialDim, Frame>& inertial_coords);
// const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
// const tnsr::ii<DataType, SpatialDim, Frame>& weyl_magnetic) {
/// @}

namespace Tags {
///
///
///
///
template <size_t SpatialDim, typename Frame, typename DataType>
struct Psi4Compute : Psi4<DataType>, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpatialRicci<SpatialDim, Frame, DataType>,
                 gr::Tags::ExtrinsicCurvature<SpatialDim, Frame, DataType>,
                 ::Tags::deriv<gr::Tags::ExtrinsicCurvature<3, Frame, DataType>,
                               tmpl::size_t<3>, Frame>,
                 gr::Tags::SpatialMetric<SpatialDim, Frame, DataType>,
                 gr::Tags::InverseSpatialMetric<SpatialDim, Frame, DataType>,
                 domain::Tags::Coordinates<SpatialDim, Frame>>;

  using return_type = Scalar<DataType>;
  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<DataType>*>,
                           const tnsr::ii<DataType, SpatialDim, Frame>&,
                           const tnsr::ii<DataType, SpatialDim, Frame>&,
                           const tnsr::ijj<DataType, SpatialDim, Frame>&,
                           const tnsr::ii<DataType, SpatialDim, Frame>&,
                           const tnsr::II<DataType, SpatialDim, Frame>&,
                           const tnsr::I<DataType, SpatialDim, Frame>&)>(
          &psi_4<SpatialDim, Frame, DataType>);
  using base = Psi4<DataType>;
};
}  // namespace Tags
}  // namespace gr
