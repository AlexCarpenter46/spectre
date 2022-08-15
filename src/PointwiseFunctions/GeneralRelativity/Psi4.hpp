// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddfef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
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
void Psi_4(const gsl::not_null<Scalar<DataType>*> psi_4_result,
           const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
           const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
           const tnsr::i<Frame>& r_hat,
           const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
           const tnsr::ii<DataType, SpatialDim, Frame>& weyl_magnetic);

template <size_t SpatialDim, typenmae Frame, typename DataType>
Scalar<DataType> Psi_4(
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::i<Frame>& r_hat,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_magnetic);
/// @}

namespace Tags {
///
///
///
///
template <size_t SpatialDim, typename Frame, typename DataType>
struct Psi4Compute : Psi4<DataType>, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<SpatialDim, Frame, DataType>,
                 StrahlkorperTags::Rhat<Frame>,
                 gr::Tags::WeylElectric<SpatialDim, Frame, DataType>,
                 gr::Tags::WelyMagnetic<DataType, 3, Frame>>;

  using return_type = Scalar<DataType>;
  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<Scalar<DataType>*>,
      const tnsr::ii<DataType, SpatialDim, Frame>&,
      const tnsr::II<DataType, SpatialDim, Frame>&, const tnsr::i<Frame>&,
      const tnsr::ii<DataType, SpatialDim, Frame>&,
      const tnsr::ii<DataType, SpatialDim, Frame>&)>(
      &Psi_4<SpatialDim, Frame, DataType>);
  using base = Psi4<DataType>;
};
}  // namespace Tags
}  // namespace gr
