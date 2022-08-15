// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/Psi4.hpp"

#include <cstddef>

#include "DataStructures/LeviCivitaIterator.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "PointwiseFunctions/GeneralRelativity/ProjectionOperators.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace gr {
template <size_t SpatialDim, typename Frame, typename DataType>
void Psi_4(const gsl::not_null<Scalar<DataType>*> psi_4_result,
           const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
           const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
           const tnsr::i<Frame>& r_hat,
           const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
           const tnsr::ii<DataType, SpatialDim, Frame>& weyl_magnetic) {
  // let's make this section for the projection tensor
  // so we'll be using the Projection Operators that are already coded up
  // cause we don't need to reinvent the wheel here.
  // They are in ProjectionOperators.hpp
  projection_tensor = transverse_projection_tensor(spatial_metric, r_hat);
  raised_r_hat = raise_or_lower_first_index(r_hat, inverse_spatial_metric);
  inverse_projection_tensor =
      transverse_projection_operator(inverse_spatial_metric, raised_r_hat);

  // projection portion of U^8+-
  tnsr::IJkk<DataType, SpatialDim, Frame> projection_product

      // weyl portion of U^8+- this is the real portion
      auto weyl_product =
  // let's make this section for the U^8+- tensor
  //
  // inside here however, we'll have to have sections because
  // this quantity is definitely a bit more complicated than the others
  // I've done so far. For instance, I think we'll to compute the
  //(P(^a_iP^b)_j - P^abP_ij) on it's own then compute the
  //(E_ab -+ e_i^kl*n_d*B_kl)
}

template <size_t SpatialDim, typenmae Frame, typename DataType>
Scalar<DataType> Psi_4(
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::i<Frame>& r_hat,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_magnetic) {}
}  // namespace gr
