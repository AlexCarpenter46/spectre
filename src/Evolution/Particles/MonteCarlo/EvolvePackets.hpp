// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"

/// Items related to the evolution of particles
/// Items related to Monte-Carlo radiation transport
namespace Particles::MonteCarlo {

struct Packet;

namespace detail {

// Time derivative of the spatial component of the momentum one-form on a null
// geodesic
void time_derivative_momentum_geodesic(
    gsl::not_null<std::array<double, 3>*> dt_momentum, const Packet& packet,
    const Scalar<DataVector>& lapse,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_lapse,
    const tnsr::iJ<DataVector, 3, Frame::Inertial>& d_shift,
    const tnsr::iJJ<DataVector, 3, Frame::Inertial>& d_inv_spatial_metric);

// Advances a single packet by time time_step along a geodesic
void evolve_single_packet_on_geodesic(
    gsl::not_null<Packet*> packet, const double& final_time,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_lapse,
    const tnsr::iJ<DataVector, 3, Frame::Inertial>& d_shift,
    const tnsr::iJJ<DataVector, 3, Frame::Inertial>& d_inv_spatial_metric,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inv_spatial_metric,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>& mesh_velocity,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian);

// Functions to be implemented to complete implementation of Monte-Carlo
// time step
double compute_fluid_frame_energy(const Packet& /*packet*/);
void compute_opacities(gsl::not_null<double*> absorption_opacity,
                       gsl::not_null<double*> scattering_opacity,
                       const double& /*fluid_frame_energy*/);
void scatter_packet(gsl::not_null<Packet*> /*packet*/);
void diffuse_packet(gsl::not_null<Packet*> /*packet*/,
                    const double& /*time_step*/);
}  // namespace detail

/// Evolve all packets in the provided std::vector (approximately) to the
/// end of the current time step.
void evolve_packets(
    gsl::not_null<std::vector<Packet>*> packets, const double& time_step,
    const Mesh<3>& mesh,
    const tnsr::I<DataVector, 3, Frame::ElementLogical>& mesh_coordinates,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_lapse,
    const tnsr::iJ<DataVector, 3, Frame::Inertial>& d_shift,
    const tnsr::iJJ<DataVector, 3, Frame::Inertial>& d_inv_spatial_metric,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inv_spatial_metric,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>& mesh_velocity,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian);

}  // namespace Particles::MonteCarlo