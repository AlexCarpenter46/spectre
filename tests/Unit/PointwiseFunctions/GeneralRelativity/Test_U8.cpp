// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <iostream>  // TODO: remove
#include <random>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Expressions/Evaluate.hpp"
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/PointwiseFunctions/GeneralRelativity/TestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/Psi4.hpp"
#include "PointwiseFunctions/GeneralRelativity/TagsDeclarations.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylPropagating.hpp"

namespace {
template <size_t SpatialDim, typename RealDataType, typename ComplexDataType>
void test_compute_item_in_databox(
    const RealDataType& used_for_size_real,
    const ComplexDataType& /*used_for_size_complex*/) {
  TestHelpers::db::test_compute_tag<gr::Tags::Psi4Compute<
      ComplexDataType, SpatialDim, Frame::Inertial, RealDataType>>("Psi4");

  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> distribution(0.1, 3.0);
  const auto nn_generator = make_not_null(&generator);
  const auto nn_distribution = make_not_null(&distribution);

  const auto spatial_ricci = make_with_random_values<
      tnsr::ii<RealDataType, SpatialDim, Frame::Inertial>>(
      nn_generator, nn_distribution, used_for_size_real);
  const auto extrinsic_curvature = make_with_random_values<
      tnsr::ii<RealDataType, SpatialDim, Frame::Inertial>>(
      nn_generator, nn_distribution, used_for_size_real);
  const auto cov_deriv_extrinsic_curvature = make_with_random_values<
      tnsr::ijj<RealDataType, SpatialDim, Frame::Inertial>>(
      nn_generator, nn_distribution, used_for_size_real);
  const auto spatial_metric =
      TestHelpers::gr::random_spatial_metric<SpatialDim, RealDataType,
                                             Frame::Inertial>(
          nn_generator, used_for_size_real);
  const auto inv_spatial_metric =
      determinant_and_inverse(spatial_metric).second;
  const auto sprt_det_spatial_metric =
      sqrt(determinant_and_inverse(spatial_metric).first);
  const auto inertial_coords = make_with_random_values<
      tnsr::I<RealDataType, SpatialDim, Frame::Inertial>>(
      nn_generator, nn_distribution, used_for_size_real);
  const auto E_ij = gr::weyl_electric(spatial_ricci, extrinsic_curvature,
                                      inverse_spatial_metric);
  const auto B_ij = gr::weyl_magnetic(cov_deriv_extrinsic_curvature,
                                      spatial_metric, sqrt_det_spatial_metric);

  const auto magnitude_cartesian =
      Scalar<RealDataType>(magnitude(inertial_coords, spatial_metric));
  auto r_hat = make_with_value<tnsr::I<RealDataType, SpatialDim, Frame>>(
      get<0, 0>(spatial_metric), 0.0);
  for (size_t j = 0; j < SpatialDim; j++) {
    // for (size_t k = 0; k < SpatialDim; k++) {
    if constexpr (is_derived_of_vector_impl_v<RealDataType>) {
      for (size_t i = 0; i < magnitude_cartesian.size(); i++) {
        if (magnitude_cartesian.get()[i] == 0.0) {
          // remove for redundancy
          r_hat.get(j)[i] = 0.0;
        } else {
          r_hat.get(j)[i] =
              inertial_coords.get(j)[i] / magnitude_cartesian.get()[i];
        }
      }
    } else {
      if (magnitude_cartesian.get() == 0.0) {
        r_hat.get(j) = 0.0;
      } else {
        r_hat.get(j) = inertial_coords.get(j) / magnitude_cartesian.get();
      }
    }
    //}
  }
  // std::cout << magnitude(r_hat, spatial_metric) << std::endl;
  const auto lower_r_hat =
      tenex::evaluate<ti::i>(r_hat(ti::J) * spatial_metric(ti::i, ti::j));
  std::cout << magnitude(lower_r_hat, inverse_spatial_metric) << std::endl;
  const auto projection_tensor =
      transverse_projection_operator(spatial_metric, lower_r_hat);
  const auto inverse_projection_tensor =
      transverse_projection_operator(inverse_spatial_metric, r_hat);
  // documentation says a and b but here since they're all spatial, we'll do
  // a = k and b = l
  const auto projection_up_lo = tenex::evaluate<ti::K, ti::i>(
      projection_tensor(ti::i, ti::j) * inverse_spatial_metric(ti::K, ti::J));
  const auto expected_U8 = gr::weyl_propagating(
      spatial_ricci, extrinsic_curvature, inverse_spatial_metric,
      cov_deriv_extrinsic_curvature, r_hat, inverse_projection_tensor,
      projection_tensor, projection_up_lo, 1.0);

  const auto new_U8 = E_ij;
  for (size_t i = 0;)
