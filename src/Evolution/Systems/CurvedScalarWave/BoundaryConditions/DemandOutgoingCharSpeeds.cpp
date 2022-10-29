// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Characteristics.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeString.hpp"

namespace CurvedScalarWave::BoundaryConditions {

template <size_t Dim>
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
DemandOutgoingCharSpeeds<Dim>::get_clone() const {
  return std::make_unique<DemandOutgoingCharSpeeds>(*this);
}

template <size_t Dim>
void DemandOutgoingCharSpeeds<Dim>::pup(PUP::er& p) {
  BoundaryCondition<Dim>::pup(p);
}

template <size_t Dim>
DemandOutgoingCharSpeeds<Dim>::DemandOutgoingCharSpeeds(
    CkMigrateMessage* const msg)
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::optional<std::string>
DemandOutgoingCharSpeeds<Dim>::dg_demand_outgoing_char_speeds(
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, Dim>& normal_covector,
    const tnsr::I<DataVector, Dim>& /*normal_vector*/,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim>& shift) const {
  tnsr::a<DataVector, 3, Frame::Inertial> char_speeds{lapse.size()};
  characteristic_speeds(make_not_null(&char_speeds), gamma1, lapse, shift,
                        normal_covector);

  if (face_mesh_velocity.has_value()) {
    const auto face_speed = dot_product(normal_covector, *face_mesh_velocity);
    for (auto& char_speed : char_speeds) {
      char_speed -= get(face_speed);
    }
  }
  for (size_t i = 0; i < char_speeds.size(); ++i) {
    if (min(char_speeds[i]) < 0.) {
      return MakeString{}
             << "Detected negative characteristic speed at boundary with "
                "outgoing char speeds boundary conditions specified. The "
                "speed is "
             << min(char_speeds[i]) << " for index " << i
             << ". To see which characteristic field this corresponds to, "
                "check the function `characteristic_speeds` in "
                "Evolution/Systems/CurvedScalarWave/Characteristics.hpp.";
    }
  }
  return std::nullopt;  // LCOV_EXCL_LINE
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID DemandOutgoingCharSpeeds<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) \
  template class DemandOutgoingCharSpeeds<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace CurvedScalarWave::BoundaryConditions
