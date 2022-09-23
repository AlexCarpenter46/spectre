// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/Psi4.hpp"

#include <cstddef>
#include <iostream>

#include "DataStructures/LeviCivitaIterator.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Expressions/Evaluate.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/Tags.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/ProjectionOperators.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylPropagating.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace gr {
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
    const tnsr::I<RealDataType, SpatialDim, Frame>& inertial_coords) {
  const auto magnitude_cartesian =
      Scalar<RealDataType>(magnitude(inertial_coords, spatial_metric));
  const auto r_hat = tenex::evaluate<ti::j>(inertial_coords(ti::I) *
                                            spatial_metric(ti::i, ti::j) /
                                            magnitude_cartesian());
  const auto projection_tensor =
      transverse_projection_operator(spatial_metric, r_hat);
  const auto raised_r_hat = tenex::evaluate<ti::I>(
      r_hat(ti::j) * inverse_spatial_metric(ti::I, ti::J));
  const auto inverse_projection_tensor =
      transverse_projection_operator(inverse_spatial_metric, raised_r_hat);
  // documentation says a and b but here since they're all spatial, we'll do
  // a = k and b = l
  const auto projection_up_lo = tenex::evaluate<ti::K, ti::i>(
      projection_tensor(ti::i, ti::j) * inverse_spatial_metric(ti::K, ti::J));

  const auto u8_plus = gr::weyl_propagating(
      spatial_ricci, extrinsic_curvature, inverse_spatial_metric,
      cov_deriv_extrinsic_curvature, raised_r_hat, inverse_projection_tensor,
      projection_tensor, projection_up_lo, 1);
  auto x_coord = make_with_value<tnsr::I<RealDataType, SpatialDim, Frame>>(
      get<0, 0>(inverse_spatial_metric), 0.0);
  x_coord.get(0) = 1.0;
  const auto magnitude_x =
      Scalar<RealDataType>(magnitude(x_coord, spatial_metric));
  // std::cout << "magnitude x: " << magnitude_x << std::endl;
  // std::cout << "x_coords: " << x_coord << std::endl;
  const auto x_hat = tenex::evaluate<ti::I>(x_coord(ti::I) / magnitude_x());
  // const auto x_hat = tenex::evaluate<ti::I>(x_coord(ti::I));
  auto y_coord = make_with_value<tnsr::I<RealDataType, SpatialDim, Frame>>(
      get<0, 0>(inverse_spatial_metric), 0.0);
  y_coord.get(1) = 1.0;
  const auto magnitude_y =
      Scalar<RealDataType>(magnitude(y_coord, spatial_metric));
  // std::cout << "magnitude y: " << magnitude_y << std::endl;
  const std::complex<double> i = std::complex<double>(0.0, 1.0);
  tnsr::I<ComplexDataType, SpatialDim, Frame> y_hat{};
  // std::complex<double> + DataVector -> ComplexDataVector
  tenex::evaluate<ti::I>(make_not_null(&y_hat),
                         i * y_coord(ti::I) / magnitude_y());
  // tenex::evaluate<ti::I>(make_not_null(&y_hat),
  //                        i * y_coord(ti::I));

  // todo rad 2 stuff
  const auto m_bar = tenex::evaluate<ti::I>((x_hat(ti::I) + y_hat(ti::I)));

  // making the real portion of psi4 only for now.
  tenex::evaluate(psi_4_result,
                  0.5 * u8_plus(ti::i, ti::j) * m_bar(ti::I) * m_bar(ti::J));
}

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
    const tnsr::I<RealDataType, SpatialDim, Frame>& inertial_coords) {
  // const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
  // const tnsr::ii<DataType, SpatialDim, Frame>& weyl_magnetic) {
  auto psi_4_result = make_with_value<Scalar<ComplexDataType>>(
      get<0, 0>(inverse_spatial_metric), 0.0);
  psi_4(make_not_null(&psi_4_result), spatial_ricci, extrinsic_curvature,
        cov_deriv_extrinsic_curvature, spatial_metric, inverse_spatial_metric,
        inertial_coords);
  return psi_4_result;
}
}  // namespace gr

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define RDTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)
#define CDTYPE(data) BOOST_PP_TUPLE_ELEM(2, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(3, data)

#define INSTANTIATE(_, data)                                                 \
  template Scalar<CDTYPE(data)> gr::psi_4(                                   \
      const tnsr::ii<RDTYPE(data), DIM(data), FRAME(data)>& spatial_ricci,   \
      const tnsr::ii<RDTYPE(data), DIM(data), FRAME(data)>&                  \
          extrinsic_curvature,                                               \
      const tnsr::ijj<RDTYPE(data), DIM(data), FRAME(data)>&                 \
          cov_deriv_extrinsic_curvature,                                     \
      const tnsr::ii<RDTYPE(data), DIM(data), FRAME(data)>& spatial_metric,  \
      const tnsr::II<RDTYPE(data), DIM(data), FRAME(data)>&                  \
          inverse_spatial_metric,                                            \
      const tnsr::I<RDTYPE(data), DIM(data), FRAME(data)>& inertial_coords); \
  template void gr::psi_4(                                                   \
      const gsl::not_null<Scalar<CDTYPE(data)>*> psi_4_result,               \
      const tnsr::ii<RDTYPE(data), DIM(data), FRAME(data)>& spatial_ricci,   \
      const tnsr::ii<RDTYPE(data), DIM(data), FRAME(data)>&                  \
          extrinsic_curvature,                                               \
      const tnsr::ijj<RDTYPE(data), DIM(data), FRAME(data)>&                 \
          cov_deriv_extrinsic_curvature,                                     \
      const tnsr::ii<RDTYPE(data), DIM(data), FRAME(data)>& spatial_metric,  \
      const tnsr::II<RDTYPE(data), DIM(data), FRAME(data)>&                  \
          inverse_spatial_metric,                                            \
      const tnsr::I<RDTYPE(data), DIM(data), FRAME(data)>& inertial_coords);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (double),
                        (std::complex<double>), (Frame::Grid, Frame::Inertial))
GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (DataVector),
                        (ComplexDataVector), (Frame::Grid, Frame::Inertial))

#undef DIM
#undef RDTYPE
#undef CDTYPE
#undef FRAME
#undef INSTANTIATE
