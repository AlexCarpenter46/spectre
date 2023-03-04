// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/FixedHashMap.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/MaxNumberOfNeighbors.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/OptionalHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace fd {
/// @{
/*!
 * \brief Computes a high-order boundary correction $G$ at the FD interface.
 *
 * The correction to the second-order boundary correction is given by
 *
 * \f{align*}{
 *  G=G^{(2)}-G^{(4)}+G^{(6)}-G^{(8)}+G^{(10)},
 * \f}
 *
 * where
 *
 *\f{align*}{
 * G^{(4)}_{j+1/2}&=\frac{1}{6}\left(G_j -2 G^{(2)} +
 *                         G_{j+1}\right), \\
 * G^{(6)}_{j+1/2}&=\frac{1}{180}\left(G_{j-1} - 9 G_j + 16 G^{(2)}
 *                         -9 G_{j+1} + G_{j+2}\right), \\
 * G^{(8)}_{j+1/2}&=\frac{1}{2100}\left(G_{j-2} - \frac{25}{3} G_{j-1}
 *                         + 50 G_j - \frac{256}{3} G^{(2)} + 50 G_{j+1}
 *                         - \frac{25}{3} G_{j+2} +G_{j+3}\right), \\
 * G^{(10)}_{j+1/2}&=\frac{1}{17640}
 *                         \left(G_{j-3} - \frac{49}{5} G_{j-2}
 *                     + 49 G_{j-1} - 245 G_j + \frac{2048}{5} G^{(2)}\right.
 *                         \nonumber \\
 *                       &\left.- 245 G_{j+1}+ 49 G_{j+2} - \frac{49}{5} G_{j+3}
 *                         + G_{j+4}\right),
 * \f}
 *
 * where
 *
 * \f{align*}{
 *  G_{j} &= F^i_j n_i^{j+1/2}, \\
 *  G_{j\pm1} &= F^i_{j\pm1} n_i^{j+1/2}, \\
 *  G_{j\pm2} &= F^i_{j\pm2} n_i^{j+1/2}, \\
 *  G_{j\pm3} &= F^i_{j\pm3} n_i^{j+1/2}, \\
 *  G_{j\pm4} &= F^i_{j\pm4} n_i^{j+1/2}.
 * \f}
 *
 * This is a generalization of the correction presented in \cite CHEN2016604.
 *
 * This high-order flux can be fed into a flux limiter, e.g. to guarantee
 * positivity.
 *
 * \note This implementation should be profiled and optimized.
 *
 * \warning This documentation is for the general case. In the restricted
 * Cartesian case we use the cell-centered flux as opposed to `G^{(4)}`, which
 * differs by a minus sign. This amounts to a minus sign change in front of the
 * $G^{(k)}$ terms in computing $G$ for $k>2$, and also a sign change in front
 * of $G^{(2)}$ in all $G^{(k)}$ for $k>2$.
 */
template <size_t CorrectionOrder, size_t Dim, typename... EvolvedVarsTags>
void cartesian_high_order_fluxes_using_nodes(
    const gsl::not_null<
        std::array<Variables<tmpl::list<EvolvedVarsTags...>>, Dim>*>
        high_order_boundary_corrections_in_logical_direction,

    const std::array<Variables<tmpl::list<EvolvedVarsTags...>>, Dim>&
        second_order_boundary_corrections_in_logical_direction,
    const Variables<tmpl::list<
        ::Tags::Flux<EvolvedVarsTags, tmpl::size_t<Dim>, Frame::Inertial>...>>&
        cell_centered_inertial_flux,
    const DirectionMap<
        Dim, Variables<tmpl::list<::Tags::Flux<
                 EvolvedVarsTags, tmpl::size_t<Dim>, Frame::Inertial>...>>>&
        ghost_cell_inertial_flux,
    const Mesh<Dim>& subcell_mesh, const size_t number_of_ghost_cells) {
  static constexpr size_t max_correction_order = 10;
  static constexpr size_t correction_width = CorrectionOrder / 2 - 1;
  static_assert(CorrectionOrder <= max_correction_order);
  for (size_t dim = 0; dim < Dim; ++dim) {
    gsl::at(*high_order_boundary_corrections_in_logical_direction, dim)
        .initialize(
            gsl::at(second_order_boundary_corrections_in_logical_direction, dim)
                .number_of_grid_points());
  }

  const auto impl = [&cell_centered_inertial_flux, &ghost_cell_inertial_flux,
                     &high_order_boundary_corrections_in_logical_direction,
                     number_of_ghost_cells,
                     &second_order_boundary_corrections_in_logical_direction,
                     &subcell_mesh](auto tag_v, auto dim_v) {
    using tag = decltype(tag_v);
    constexpr size_t dim = decltype(dim_v)::value;

    auto& high_order_var_correction =
        get<tag>((*high_order_boundary_corrections_in_logical_direction)[dim]);
    const auto& second_order_var_correction =
        get<tag>(second_order_boundary_corrections_in_logical_direction[dim]);
    const auto& cell_centered_flux =
        get<::Tags::Flux<tag, tmpl::size_t<Dim>, Frame::Inertial>>(
            cell_centered_inertial_flux);
    const auto& lower_neighbor_cell_centered_flux =
        get<::Tags::Flux<tag, tmpl::size_t<Dim>, Frame::Inertial>>(
            ghost_cell_inertial_flux.at(Direction<Dim>{dim, Side::Lower}));
    const auto& upper_neighbor_cell_centered_flux =
        get<::Tags::Flux<tag, tmpl::size_t<Dim>, Frame::Inertial>>(
            ghost_cell_inertial_flux.at(Direction<Dim>{dim, Side::Upper}));
    using FluxTensor = std::decay_t<decltype(cell_centered_flux)>;
    const auto& subcell_extents = subcell_mesh.extents();
    auto subcell_face_extents = subcell_extents;
    ++subcell_face_extents[dim];
    auto neighbor_extents = subcell_extents;
    neighbor_extents[dim] = number_of_ghost_cells;
    const size_t number_of_components = second_order_var_correction.size();
    for (size_t storage_index = 0; storage_index < number_of_components;
         ++storage_index) {
      const auto flux_multi_index = prepend(
          second_order_var_correction.get_tensor_index(storage_index), dim);
      const size_t flux_storage_index =
          FluxTensor::get_storage_index(flux_multi_index);
      // Loop over each face
      for (size_t k = 0; k < (Dim == 3 ? subcell_face_extents[2] : 1); ++k) {
        for (size_t j = 0; j < (Dim >= 2 ? subcell_face_extents[1] : 1); ++j) {
          for (size_t i = 0; i < subcell_face_extents[0]; ++i) {
            const Index<Dim> face_index = [i, j, k]() -> Index<Dim> {
              if constexpr (Dim == 3) {
                return Index<Dim>{i, j, k};
              } else if constexpr (Dim == 2) {
                (void)k;
                return Index<Dim>{i, j};
              } else {
                (void)k, (void)j;
                return Index<Dim>{i};
              }
            }();
            const size_t face_storage_index =
                collapsed_index(face_index, subcell_face_extents);
            Index<Dim> neighbor_index{};
            for (size_t l = 0; l < Dim; ++l) {
              if (l != dim) {
                neighbor_index[l] = face_index[l];
              }
            }

            double& correction =
                high_order_var_correction[storage_index][face_storage_index] =
                    0.0;

            std::array<double, CorrectionOrder - 2>
                cell_centered_fluxes_for_stencil{};
            // fill if we have to retrieve from lower neighbor
            size_t stencil_index = 0;
            for (int grid_index = static_cast<int>(face_index[dim]) -
                                  static_cast<int>(correction_width);
                 grid_index < static_cast<int>(face_index[dim]) +
                                  static_cast<int>(correction_width);
                 ++grid_index, ++stencil_index) {
              if (grid_index < 0) {
                neighbor_index[dim] = static_cast<size_t>(
                    static_cast<int>(number_of_ghost_cells) + grid_index);
                gsl::at(cell_centered_fluxes_for_stencil, stencil_index) =
                    lower_neighbor_cell_centered_flux[flux_storage_index]
                                                     [collapsed_index(
                                                         neighbor_index,
                                                         neighbor_extents)];
              } else if (grid_index >= static_cast<int>(subcell_extents[dim])) {
                neighbor_index[dim] = static_cast<size_t>(
                    grid_index - static_cast<int>(subcell_extents[dim]));
                gsl::at(cell_centered_fluxes_for_stencil, stencil_index) =
                    upper_neighbor_cell_centered_flux[flux_storage_index]
                                                     [collapsed_index(
                                                         neighbor_index,
                                                         neighbor_extents)];
              } else {
                Index<Dim> volume_index = face_index;
                volume_index[dim] = static_cast<size_t>(grid_index);
                gsl::at(cell_centered_fluxes_for_stencil, stencil_index) =
                    cell_centered_flux[flux_storage_index][collapsed_index(
                        volume_index, subcell_extents)];
              }
            }

            // Only do 4th-order correction for now...
            if constexpr (CorrectionOrder >= 10) {
              correction -=
                  5.6689342403628117913e-5 *
                  (gsl::at(cell_centered_fluxes_for_stencil,
                           correction_width - 4) +
                   gsl::at(cell_centered_fluxes_for_stencil,
                           correction_width + 3) -
                   9.8 * (gsl::at(cell_centered_fluxes_for_stencil,
                                  correction_width - 3) +
                          gsl::at(cell_centered_fluxes_for_stencil,
                                  correction_width + 2)) +
                   49.0 * (gsl::at(cell_centered_fluxes_for_stencil,
                                   correction_width - 2) +
                           gsl::at(cell_centered_fluxes_for_stencil,
                                   correction_width + 1)) -
                   245.0 * (gsl::at(cell_centered_fluxes_for_stencil,
                                    correction_width - 1) +
                            gsl::at(cell_centered_fluxes_for_stencil,
                                    correction_width)) -
                   409.6 * second_order_var_correction[storage_index]
                                                      [face_storage_index]);
            }
            if constexpr (CorrectionOrder >= 8) {
              correction +=
                  4.7619047619047619047e-4 *
                  (gsl::at(cell_centered_fluxes_for_stencil,
                           correction_width - 3) +
                   gsl::at(cell_centered_fluxes_for_stencil,
                           correction_width + 2) -
                   8.3333333333333333333 *
                       (gsl::at(cell_centered_fluxes_for_stencil,
                                correction_width - 2) +
                        gsl::at(cell_centered_fluxes_for_stencil,
                                correction_width + 1)) +
                   50.0 * (gsl::at(cell_centered_fluxes_for_stencil,
                                   correction_width - 1) +
                           gsl::at(cell_centered_fluxes_for_stencil,
                                   correction_width)) +
                   85.333333333333333333 *
                       second_order_var_correction[storage_index]
                                                  [face_storage_index]);
            }
            if constexpr (CorrectionOrder >= 6) {
              correction -=
                  5.5555555555555555555e-3 *
                  (gsl::at(cell_centered_fluxes_for_stencil,
                           correction_width - 2) +
                   gsl::at(cell_centered_fluxes_for_stencil,
                           correction_width + 1) -
                   9.0 * (gsl::at(cell_centered_fluxes_for_stencil,
                                  correction_width - 1) +
                          gsl::at(cell_centered_fluxes_for_stencil,
                                  correction_width)) -
                   16.0 * second_order_var_correction[storage_index]
                                                     [face_storage_index]);
            }
            if constexpr (CorrectionOrder >= 4) {
              correction +=
                  0.166666666666666666 *
                  (gsl::at(cell_centered_fluxes_for_stencil,
                           correction_width - 1) +
                   gsl::at(cell_centered_fluxes_for_stencil, correction_width) +
                   2.0 * second_order_var_correction[storage_index]
                                                    [face_storage_index]);
            }

            // Add second-order correction last
            correction +=
                second_order_var_correction[storage_index][face_storage_index];
          }
        }
      }
    }
  };

  EXPAND_PACK_LEFT_TO_RIGHT(
      impl(EvolvedVarsTags{}, std::integral_constant<size_t, 0>{}));
  if constexpr (Dim > 1) {
    EXPAND_PACK_LEFT_TO_RIGHT(
        impl(EvolvedVarsTags{}, std::integral_constant<size_t, 1>{}));
    if constexpr (Dim > 2) {
      EXPAND_PACK_LEFT_TO_RIGHT(
          impl(EvolvedVarsTags{}, std::integral_constant<size_t, 2>{}));
    }
  }
}

template <size_t Dim, typename... EvolvedVarsTags>
void cartesian_high_order_fluxes_using_nodes(
    const gsl::not_null<
        std::array<Variables<tmpl::list<EvolvedVarsTags...>>, Dim>*>
        high_order_boundary_corrections_in_logical_direction,

    const std::array<Variables<tmpl::list<EvolvedVarsTags...>>, Dim>&
        second_order_boundary_corrections_in_logical_direction,
    const Variables<tmpl::list<
        ::Tags::Flux<EvolvedVarsTags, tmpl::size_t<Dim>, Frame::Inertial>...>>&
        cell_centered_inertial_flux,
    const DirectionMap<
        Dim, Variables<tmpl::list<::Tags::Flux<
                 EvolvedVarsTags, tmpl::size_t<Dim>, Frame::Inertial>...>>>&
        ghost_cell_inertial_flux,
    const Mesh<Dim>& subcell_mesh, const size_t number_of_ghost_cells,
    const size_t correction_order) {
  switch (correction_order) {
    case 2:
      cartesian_high_order_fluxes_using_nodes<2>(
          high_order_boundary_corrections_in_logical_direction,
          second_order_boundary_corrections_in_logical_direction,
          cell_centered_inertial_flux, ghost_cell_inertial_flux, subcell_mesh,
          number_of_ghost_cells);
      break;
    case 4:
      cartesian_high_order_fluxes_using_nodes<4>(
          high_order_boundary_corrections_in_logical_direction,
          second_order_boundary_corrections_in_logical_direction,
          cell_centered_inertial_flux, ghost_cell_inertial_flux, subcell_mesh,
          number_of_ghost_cells);
      break;
    case 6:
      cartesian_high_order_fluxes_using_nodes<6>(
          high_order_boundary_corrections_in_logical_direction,
          second_order_boundary_corrections_in_logical_direction,
          cell_centered_inertial_flux, ghost_cell_inertial_flux, subcell_mesh,
          number_of_ghost_cells);
      break;
    case 8:
      cartesian_high_order_fluxes_using_nodes<8>(
          high_order_boundary_corrections_in_logical_direction,
          second_order_boundary_corrections_in_logical_direction,
          cell_centered_inertial_flux, ghost_cell_inertial_flux, subcell_mesh,
          number_of_ghost_cells);
      break;
    case 10:
      cartesian_high_order_fluxes_using_nodes<10>(
          high_order_boundary_corrections_in_logical_direction,
          second_order_boundary_corrections_in_logical_direction,
          cell_centered_inertial_flux, ghost_cell_inertial_flux, subcell_mesh,
          number_of_ghost_cells);
      break;
    default:
      ERROR("Unsupported correction order " << correction_order
                                            << ". Only know 2,4,6,8,10");
  };
}
/// @}

/*!
 * \brief Fill the `flux_neighbor_data` with pointers into the
 * `reconstruction_neighbor_data`.
 *
 * The `reconstruction_neighbor_data` is stored in the tag
 * `evolution::dg::subcell::Tags::NeighborDataForReconstruction`, and the
 * `ghost_zone_size` should come from the FD reconstructor.
 */
template <size_t Dim, typename FluxesTags>
void set_cartesian_neighbor_cell_centered_fluxes(
    const gsl::not_null<DirectionMap<Dim, Variables<FluxesTags>>*>
        flux_neighbor_data,
    const FixedHashMap<maximum_number_of_neighbors(Dim),
                       std::pair<Direction<Dim>, ElementId<Dim>>, DataVector,
                       boost::hash<std::pair<Direction<Dim>, ElementId<Dim>>>>&
        reconstruction_neighbor_data,
    const Mesh<Dim>& subcell_mesh, const size_t ghost_zone_size,
    const size_t number_of_rdmp_values_in_ghost_data) {
  for (const auto& [direction_id, ghost_data] : reconstruction_neighbor_data) {
    const size_t neighbor_flux_size =
        subcell_mesh.number_of_grid_points() /
        subcell_mesh.extents(direction_id.first.dimension()) * ghost_zone_size *
        Variables<FluxesTags>::number_of_independent_components;
    (*flux_neighbor_data)[direction_id.first].set_data_ref(
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
        const_cast<double*>(std::next(
            ghost_data.data(),
            static_cast<std::ptrdiff_t>(ghost_data.size() -
                                        number_of_rdmp_values_in_ghost_data -
                                        neighbor_flux_size))),
        neighbor_flux_size);
  }
}

/*!
 * \brief Computes the high-order Cartesian flux corrections if necessary.
 *
 * The `cell_centered_fluxes` is stored in the tag
 * `evolution::dg::subcell::Tags::CellCenteredFlux`, `fd_derivative_order` is
 * from `evolution::dg::subcell::Tags::SubcellOptions`
 * (`.finite_difference_derivative_order()`), the `reconstruction_neighbor_data`
 * is stored in the tag
 * `evolution::dg::subcell::Tags::NeighborDataForReconstruction`, the
 * `ghost_zone_size` should come from the FD reconstructor.
 *
 * By default we assume no RDMP data is in the `ghost_data` buffer. In the
 * future we will want to update how we store the data in order to eliminate
 * more memory allocations and copies, in which case that value will be
 * non-zero.
 *
 * \note `high_order_corrections` must either not have a value or have all
 * elements be of the same size as
 * `second_order_boundary_corrections[0].number_of_grid_points()`, where we've
 * assumed `second_order_boundary_corrections` is the same in all directions.
 */
template <size_t Dim, typename... EvolvedVarsTags,
          typename FluxesTags = tmpl::list<::Tags::Flux<
              EvolvedVarsTags, tmpl::size_t<Dim>, Frame::Inertial>...>>
void cartesian_high_order_flux_corrections(
    const gsl::not_null<std::optional<
        std::array<Variables<tmpl::list<EvolvedVarsTags...>>, Dim>>*>
        high_order_corrections,

    const std::optional<Variables<FluxesTags>>& cell_centered_fluxes,
    const std::array<Variables<tmpl::list<EvolvedVarsTags...>>, Dim>&
        second_order_boundary_corrections,
    const std::optional<size_t>& fd_derivative_order,
    const FixedHashMap<maximum_number_of_neighbors(Dim),
                       std::pair<Direction<Dim>, ElementId<Dim>>, DataVector,
                       boost::hash<std::pair<Direction<Dim>, ElementId<Dim>>>>&
        reconstruction_neighbor_data,
    const Mesh<Dim>& subcell_mesh, const size_t ghost_zone_size,
    const size_t number_of_rdmp_values_in_ghost_data = 0) {
  if (cell_centered_fluxes.has_value()) {
    ASSERT(fd_derivative_order.has_value(),
           "No finite difference order is set but we are trying to do "
           "high-order FD.");
    ASSERT(alg::all_of(
               second_order_boundary_corrections,
               [expected_size = second_order_boundary_corrections[0]
                                    .number_of_grid_points()](const auto& e) {
                 return e.number_of_grid_points() == expected_size;
               }),
           "All second-order boundary corrections must be of the same size, "
               << second_order_boundary_corrections[0].number_of_grid_points());
    if (fd_derivative_order.value_or(2) > 2) {
      if (not high_order_corrections->has_value()) {
        (*high_order_corrections) =
            make_array<Dim>(Variables<tmpl::list<EvolvedVarsTags...>>{
                second_order_boundary_corrections[0].number_of_grid_points()});
      }
      ASSERT(
          high_order_corrections->has_value() and
              alg::all_of(high_order_corrections->value(),
                          [expected_size =
                               second_order_boundary_corrections[0]
                                   .number_of_grid_points()](const auto& e) {
                            return e.number_of_grid_points() == expected_size;
                          }),
          "The high_order_corrections must all have size "
              << second_order_boundary_corrections[0].number_of_grid_points());
      DirectionMap<Dim, Variables<FluxesTags>> flux_neighbor_data{};
      set_cartesian_neighbor_cell_centered_fluxes(
          make_not_null(&flux_neighbor_data), reconstruction_neighbor_data,
          subcell_mesh, ghost_zone_size, number_of_rdmp_values_in_ghost_data);

      cartesian_high_order_fluxes_using_nodes(
          make_not_null(&(high_order_corrections->value())),
          second_order_boundary_corrections, cell_centered_fluxes.value(),
          flux_neighbor_data, subcell_mesh, ghost_zone_size,
          fd_derivative_order.value());
    }
  }
}
}  // namespace fd
