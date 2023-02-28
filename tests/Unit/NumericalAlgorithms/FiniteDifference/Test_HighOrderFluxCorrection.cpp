// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <random>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Evolution/DgSubcell/CartesianFluxDivergence.hpp"
#include "Evolution/DgSubcell/SliceData.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/FiniteDifference/HighOrderFluxCorrection.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"

namespace {
struct Scalar0 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

template <size_t Dim>
struct Vector0 : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::Inertial>;
};

template <size_t Dim, size_t CorrectionOrder>
void test() {
  const size_t points_per_dimension = CorrectionOrder + 2;
  CAPTURE(points_per_dimension);
  CAPTURE(CorrectionOrder);
  CAPTURE(Dim);
  const size_t max_degree = CorrectionOrder;
  const size_t stencil_width = CorrectionOrder + 1;
  const size_t number_of_ghost_points = (stencil_width - 1) / 2 + 1;

  using FluxTags = tmpl::list<Scalar0, Vector0<Dim>>;
  using Scalar0Flux = ::Tags::Flux<Scalar0, tmpl::size_t<Dim>, Frame::Inertial>;
  using Vector0Flux =
      ::Tags::Flux<Vector0<Dim>, tmpl::size_t<Dim>, Frame::Inertial>;
  using FluxVars =
      Variables<db::wrap_tags_in<::Tags::Flux, FluxTags, tmpl::size_t<Dim>,
                                 Frame::Inertial>>;
  using CorrectionVars = Variables<FluxTags>;

  const Mesh<Dim> mesh{points_per_dimension, Spectral::Basis::FiniteDifference,
                       Spectral::Quadrature::CellCentered};
  auto logical_coords = logical_coordinates(mesh);
  // Make the logical coordinates different in each direction
  for (size_t i = 1; i < Dim; ++i) {
    logical_coords.get(i) += 4.0 * static_cast<double>(i);
  }

  // Compute polynomial on cell centers in FD cluster of points
  const auto set_polynomial = Overloader{
      [max_degree](const gsl::not_null<FluxVars*> vars_ptr,
                   const auto& local_logical_coords) {
        (void)max_degree;
        for (size_t storage_index = 0;
             storage_index < get<Scalar0Flux>(*vars_ptr).size();
             ++storage_index) {
          get<Scalar0Flux>(*vars_ptr)[storage_index] = 0.0;
          for (size_t degree = 1; degree <= max_degree; ++degree) {
            for (size_t i = 0; i < Dim; ++i) {
              get<Scalar0Flux>(*vars_ptr)[storage_index] +=
                  pow(local_logical_coords.get(i), degree);
            }
          }
        }
        for (size_t storage_index = 0;
             storage_index < get<Vector0Flux>(*vars_ptr).size();
             ++storage_index) {
          get<Vector0Flux>(*vars_ptr)[storage_index] =
              1.0 + 0.3 * static_cast<double>(storage_index);
          for (size_t degree = 1; degree <= max_degree; ++degree) {
            for (size_t i = 0; i < Dim; ++i) {
              get<Vector0Flux>(*vars_ptr)[storage_index] +=
                  pow(local_logical_coords.get(i), degree);
            }
          }
        }
      },
      [max_degree](const gsl::not_null<CorrectionVars*> vars_ptr,
                   const auto& local_logical_coords) {
        (void)max_degree;
        for (size_t storage_index = 0;
             storage_index < get<Scalar0>(*vars_ptr).size(); ++storage_index) {
          get<Scalar0>(*vars_ptr)[storage_index] = 0.0;
          for (size_t degree = 1; degree <= max_degree; ++degree) {
            for (size_t i = 0; i < Dim; ++i) {
              get<Scalar0>(*vars_ptr)[storage_index] +=
                  pow(local_logical_coords.get(i), degree);
            }
          }
        }
        for (size_t storage_index = 0;
             storage_index < get<Vector0<Dim>>(*vars_ptr).size();
             ++storage_index) {
          get<Vector0<Dim>>(*vars_ptr)[storage_index] =
              100.0 + 11.0 * static_cast<double>(storage_index);
          for (size_t degree = 1; degree <= max_degree; ++degree) {
            for (size_t i = 0; i < Dim; ++i) {
              get<Vector0<Dim>>(*vars_ptr)[storage_index] +=
                  pow(local_logical_coords.get(i), degree);
            }
          }
        }
      }};
  const auto set_polynomial_divergence =
      [max_degree](const gsl::not_null<CorrectionVars*> d_vars_ptr,
                   const auto& local_logical_coords) {
        (void)max_degree;
        get(get<Scalar0>(*d_vars_ptr)) = 0.0;
        for (size_t i = 0; i < Dim; ++i) {
          // constant deriv is zero
          get<Vector0<Dim>>(*d_vars_ptr).get(i) = 0.0;
        }
        // Compute divergence
        for (size_t deriv_dim = 0; deriv_dim < Dim; ++deriv_dim) {
          for (size_t degree = 1; degree <= max_degree; ++degree) {
            get(get<Scalar0>(*d_vars_ptr)) +=
                degree * pow(local_logical_coords.get(deriv_dim), degree - 1);
            for (size_t i = 0; i < Dim; ++i) {
              get<Vector0<Dim>>(*d_vars_ptr).get(i) +=
                  degree * pow(local_logical_coords.get(deriv_dim), degree - 1);
            }
          }
        }
      };
  FluxVars volume_vars(mesh.number_of_grid_points());
  set_polynomial(&volume_vars, logical_coords);

  CorrectionVars expected_divergence(mesh.number_of_grid_points());
  set_polynomial_divergence(&expected_divergence, logical_coords);

  // Compute the polynomial at the cell center for the neighbor data that we
  // "received".
  //
  // We do this by computing the solution in our entire neighbor, then using
  // slice_data to get the subset of points that are needed.
  DirectionMap<Dim, FluxVars> neighbor_data{};
  for (const auto& direction : Direction<Dim>::all_directions()) {
    auto neighbor_logical_coords = logical_coords;
    neighbor_logical_coords.get(direction.dimension()) +=
        direction.sign() * 2.0;
    FluxVars neighbor_vars(mesh.number_of_grid_points(), 0.0);
    set_polynomial(&neighbor_vars, neighbor_logical_coords);

    DirectionMap<Dim, bool> directions_to_slice{};
    directions_to_slice[direction.opposite()] = true;
    const auto sliced_data = evolution::dg::subcell::slice_data(
        neighbor_vars, mesh.extents(), number_of_ghost_points,
        directions_to_slice, 0);
    CAPTURE(number_of_ghost_points);
    REQUIRE(sliced_data.size() == 1);
    REQUIRE(sliced_data.contains(direction.opposite()));
    REQUIRE(sliced_data.at(direction.opposite()).size() %
                FluxVars::number_of_independent_components ==
            0);
    neighbor_data[direction].initialize(
        sliced_data.at(direction.opposite()).size() /
        FluxVars::number_of_independent_components);
    std::copy(sliced_data.at(direction.opposite()).begin(),
              sliced_data.at(direction.opposite()).end(),
              neighbor_data[direction].data());
  }

  std::array<CorrectionVars, Dim> second_order_corrections{};
  for (size_t i = 0; i < Dim; ++i) {
    // Compare to analytic solution on the faces.
    const auto basis = make_array<Dim>(Spectral::Basis::FiniteDifference);
    auto quadrature = make_array<Dim>(Spectral::Quadrature::CellCentered);
    auto extents = make_array<Dim>(points_per_dimension);
    gsl::at(extents, i) = points_per_dimension + 1;
    gsl::at(quadrature, i) = Spectral::Quadrature::FaceCentered;
    const Mesh<Dim> face_centered_mesh{extents, basis, quadrature};
    auto face_logical_coords = logical_coordinates(face_centered_mesh);
    for (size_t j = 1; j < Dim; ++j) {
      face_logical_coords.get(j) += 4.0 * static_cast<double>(j);
    }
    gsl::at(second_order_corrections, i)
        .initialize(face_centered_mesh.number_of_grid_points());
    set_polynomial(make_not_null(&gsl::at(second_order_corrections, i)),
                   face_logical_coords);
    // We use n_i F^i in the code, so need to negate to get sign to agree.
    gsl::at(second_order_corrections, i) *= -1.0;
  }

  std::array<CorrectionVars, Dim> high_order_corrections = make_array<Dim>(
      CorrectionVars{second_order_corrections[0].number_of_grid_points()});
  fd::cartesian_high_order_fluxes_using_nodes<CorrectionOrder>(
      make_not_null(&high_order_corrections), second_order_corrections,
      volume_vars, neighbor_data, mesh, number_of_ghost_points);

  // Now compute the Cartesian derivative of the high_order_corrections to
  // verify that it is computed sufficiently accurately.
  const DataVector inv_jacobian{mesh.number_of_grid_points(), 1.0};
  CorrectionVars flux_divergence{mesh.number_of_grid_points(), 0.0};
  for (size_t d = 0; d < Dim; ++d) {
    // Note: assumes isotropic mesh
    const double one_over_delta_xi =
        -1.0 / (logical_coords.get(0)[1] - logical_coords.get(0)[0]);
    evolution::dg::subcell::add_cartesian_flux_divergence(
        make_not_null(&get(get<Scalar0>(flux_divergence))), one_over_delta_xi,
        inv_jacobian, get(get<Scalar0>(gsl::at(high_order_corrections, d))),
        mesh.extents(), d);
    for (size_t i = 0; i < Dim; ++i) {
      evolution::dg::subcell::add_cartesian_flux_divergence(
          make_not_null(&get<Vector0<Dim>>(flux_divergence).get(i)),
          one_over_delta_xi, inv_jacobian,
          get<Vector0<Dim>>(gsl::at(high_order_corrections, d)).get(i),
          mesh.extents(), d);
    }
  }

  // With high-order corrections roundoff can accumulate.
  Approx custom_approx = Approx::custom().epsilon(5.e-12);
  CHECK_ITERABLE_CUSTOM_APPROX(get<Scalar0>(flux_divergence),
                               get<Scalar0>(expected_divergence),
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<Vector0<Dim>>(flux_divergence),
                               get<Vector0<Dim>>(expected_divergence),
                               custom_approx);
}

SPECTRE_TEST_CASE("Unit.FiniteDifference.CartesianHighOrderFluxCorrection",
                  "[Unit][NumericalAlgorithms]") {
  test<1, 2>();
  test<1, 4>();
  test<1, 6>();
  test<1, 8>();
  test<1, 10>();

  test<2, 2>();
  test<2, 4>();
  test<2, 6>();
  test<2, 8>();
  test<2, 10>();

  test<3, 2>();
  test<3, 4>();
  test<3, 6>();
  test<3, 8>();
  test<3, 10>();
}
}  // namespace
