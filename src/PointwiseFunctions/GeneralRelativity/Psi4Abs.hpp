// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/CoordinateMaps/Tags.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Psi4.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"

/// \cond
namespace gsl {
template <typename>
struct not_null;
}  // namespace gsl
/// \endcond

namespace gr {

template <size_t SpatialDim, typename Frame, typename RealDataType>
void psi_4_abs(
    const gsl::not_null<Scalar<RealDataType>*> psi_4_abs_result,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& extrinsic_curvature,
    const tnsr::ijj<RealDataType, SpatialDim, Frame>&
        cov_deriv_extrinsic_curvature,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<RealDataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::I<RealDataType, SpatialDim, Frame>& inertial_coords) {
  const auto psi_4_complex = psi_4<ComplexDataVector>(
      spatial_ricci, extrinsic_curvature, cov_deriv_extrinsic_curvature,
      spatial_metric, inverse_spatial_metric, inertial_coords);
  get(*psi_4_abs_result) = real(get(psi_4_complex));
}

template <size_t SpatialDim, typename Frame, typename RealDataType>
Scalar<RealDataType> psi_4_abs(
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& extrinsic_curvature,
    const tnsr::ijj<RealDataType, SpatialDim, Frame>&
        cov_deriv_extrinsic_curvature,
    const tnsr::ii<RealDataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<RealDataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::I<RealDataType, SpatialDim, Frame>& inertial_coords) {
  auto psi_4_abs_result = make_with_value<Scalar<RealDataType>>(
      get<0, 0>(inverse_spatial_metric), 0.0);
  // for(size_t i = 0; i < get(psi_4_complex).size(); i++) {
  //     get(psi_4_abs_result)[i] = abs(get(psi_4_complex)[i]);
  // }
  psi_4_abs(make_not_null(&psi_4_abs_result), spatial_ricci,
            extrinsic_curvature, cov_deriv_extrinsic_curvature, spatial_metric,
            inverse_spatial_metric, inertial_coords);
  return psi_4_abs_result;
}

namespace Tags {
template <size_t SpatialDim, typename Frame, typename RealDataType>
struct Psi4AbsCompute : Psi4Abs<RealDataType>, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::SpatialRicci<SpatialDim, Frame, RealDataType>,
      gr::Tags::ExtrinsicCurvature<SpatialDim, Frame, RealDataType>,
      ::Tags::deriv<gr::Tags::ExtrinsicCurvature<3, Frame, RealDataType>,
                    tmpl::size_t<3>, Frame>,
      gr::Tags::SpatialMetric<SpatialDim, Frame, RealDataType>,
      gr::Tags::InverseSpatialMetric<SpatialDim, Frame, RealDataType>,
      domain::Tags::Coordinates<SpatialDim, Frame>>;

  using return_type = Scalar<RealDataType>;
  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<RealDataType>*>,
                           const tnsr::ii<RealDataType, SpatialDim, Frame>&,
                           const tnsr::ii<RealDataType, SpatialDim, Frame>&,
                           const tnsr::ijj<RealDataType, SpatialDim, Frame>&,
                           const tnsr::ii<RealDataType, SpatialDim, Frame>&,
                           const tnsr::II<RealDataType, SpatialDim, Frame>&,
                           const tnsr::I<RealDataType, SpatialDim, Frame>&)>(
          &psi_4_abs<SpatialDim, Frame, RealDataType>);
  using base = Psi4Abs<RealDataType>;
};
}  // namespace Tags
}  // namespace gr
