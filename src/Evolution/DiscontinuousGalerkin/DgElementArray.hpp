// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <unordered_set>
#include <vector>

#include "Domain/Block.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/DiagnosticInfo.hpp"
#include "Domain/Domain.hpp"
#include "Domain/MinimumGridSpacing.hpp"
#include "Domain/OptionTags.hpp"
#include "Domain/Structure/CreateInitialMesh.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Domain/WeightedElementDistribution.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/Printf.hpp"
#include "Parallel/Tags/ResourceInfo.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/CreateHasStaticMemberVariable.hpp"

namespace detail {
CREATE_HAS_STATIC_MEMBER_VARIABLE(use_z_order_distribution)
CREATE_HAS_STATIC_MEMBER_VARIABLE_V(use_z_order_distribution)
}  // namespace detail

/*!
 * \brief The parallel component responsible for managing the DG elements that
 * compose the computational domain
 *
 * This parallel component will perform the actions specified by the
 * `PhaseDepActionList`.
 *
 * The element assignment to processors is performed by
 * `domain::BlockZCurveProcDistribution` (using a Morton space-filling curve),
 * unless `static constexpr bool use_z_order_distribution = false;` is specified
 * in the `Metavariables`, in which case elements are assigned to processors via
 * round-robin assignment. In both cases, an unordered set of `size_t`s can be
 * passed to the `allocate_array` function which represents physical processors
 * to avoid placing elements on.
 */
template <class Metavariables, class PhaseDepActionList>
struct DgElementArray {
  static constexpr size_t volume_dim = Metavariables::volume_dim;

  using chare_type = Parallel::Algorithms::Array;
  using metavariables = Metavariables;
  using phase_dependent_action_list = PhaseDepActionList;
  using array_index = ElementId<volume_dim>;

  using const_global_cache_tags = tmpl::list<domain::Tags::Domain<volume_dim>>;

  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
          initialization_items,
      const std::unordered_set<size_t>& procs_to_ignore = {});

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<DgElementArray>(local_cache)
        .start_phase(next_phase);
  }
};

template <class Metavariables, class PhaseDepActionList>
void DgElementArray<Metavariables, PhaseDepActionList>::allocate_array(
    Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
    const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
        initialization_items,
    const std::unordered_set<size_t>& procs_to_ignore) {
  auto& local_cache = *Parallel::local_branch(global_cache);
  auto& dg_element_array =
      Parallel::get_parallel_component<DgElementArray>(local_cache);
  const auto& domain =
      Parallel::get<domain::Tags::Domain<volume_dim>>(local_cache);
  const auto& initial_refinement_levels =
      get<domain::Tags::InitialRefinementLevels<volume_dim>>(
          initialization_items);
  const auto& initial_extents =
      get<domain::Tags::InitialExtents<volume_dim>>(initialization_items);
  const auto& quadrature =
      get<evolution::dg::Tags::Quadrature>(initialization_items);

  bool use_z_order_distribution = true;
  if constexpr (detail::has_use_z_order_distribution_v<Metavariables>) {
    use_z_order_distribution = Metavariables::use_z_order_distribution;
  }

  const size_t number_of_procs = Parallel::number_of_procs<size_t>(local_cache);
  const size_t number_of_nodes = Parallel::number_of_nodes<size_t>(local_cache);

  // Will be used to print domain diagnostic info
  std::vector<size_t> elements_per_core(number_of_procs, 0_st);
  std::vector<size_t> elements_per_node(number_of_nodes, 0_st);
  std::vector<size_t> grid_points_per_core(number_of_procs, 0_st);
  std::vector<size_t> grid_points_per_node(number_of_nodes, 0_st);
  std::vector<double> cost_per_core(number_of_procs, 0_st);
  std::vector<double> cost_per_node(number_of_nodes, 0_st);

  if (use_z_order_distribution) {
    std::vector<std::vector<double>> cost_by_element_by_block(
        domain.blocks().size());

    for (size_t block_number = 0; block_number < domain.blocks().size();
         block_number++) {
      const auto& block = domain.blocks()[block_number];
      const auto initial_ref_levs = initial_refinement_levels[block.id()];
      const std::vector<ElementId<volume_dim>> element_ids =
          initial_element_ids_in_z_score_order(block.id(), initial_ref_levs);
      const size_t grid_points_per_element = alg::accumulate(
          initial_extents[block.id()], 1_st, std::multiplies<size_t>());

      cost_by_element_by_block[block_number].reserve(element_ids.size());

      for (const auto& element_id : element_ids) {
        // TODO : move this out of here, probably best to put it in
        // WeightedElementDistribution to keep all the logic handled there in
        // one class.
        Mesh<volume_dim> mesh = ::domain::Initialization::create_initial_mesh(
            initial_extents, element_id, quadrature);
        Element<volume_dim> element =
            ::domain::Initialization::create_initial_element(
                element_id, block, initial_refinement_levels);
        ElementMap<volume_dim, Frame::Grid> element_map{
            element_id,
            block.is_time_dependent()
                ? block.moving_mesh_logical_to_grid_map().get_clone()
                : block.stationary_map().get_to_grid_frame()};

        std::unique_ptr<::domain::CoordinateMapBase<
            Frame::Grid, Frame::Inertial, volume_dim>>
            grid_to_inertial_map;
        if (block.is_time_dependent()) {
          grid_to_inertial_map =
              block.moving_mesh_grid_to_inertial_map().get_clone();
        } else {
          grid_to_inertial_map =
              ::domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
                  ::domain::CoordinateMaps::Identity<volume_dim>{});
        }

        tnsr::I<DataVector, volume_dim, Frame::ElementLogical> logical_coords{};
        domain::Tags::LogicalCoordinates<volume_dim>::function(
            make_not_null(&logical_coords), mesh);

        tnsr::I<DataVector, volume_dim, Frame::Grid> grid_coords{};
        domain::Tags::MappedCoordinates<
            domain::Tags::ElementMap<volume_dim, Frame::Grid>,
            domain::Tags::Coordinates<volume_dim, Frame::ElementLogical>>::
            function(make_not_null(&grid_coords), element_map, logical_coords);

        double minimum_grid_spacing =
            std::numeric_limits<double>::signaling_NaN();
        domain::Tags::MinimumGridSpacingCompute<volume_dim, Frame::Grid>::
            function(make_not_null(&minimum_grid_spacing), mesh, grid_coords);

        cost_by_element_by_block[block_number].emplace_back(
            grid_points_per_element / sqrt(minimum_grid_spacing));
      }
    }

    const size_t num_of_procs_to_use = number_of_procs - procs_to_ignore.size();
    const domain::WeightedBlockZCurveProcDistribution<volume_dim>
        element_distribution{num_of_procs_to_use, cost_by_element_by_block,
                             procs_to_ignore};

    std::vector<size_t> grid_points_by_element{};

    for (size_t block_number = 0; block_number < domain.blocks().size();
         block_number++) {
      const auto& block = domain.blocks()[block_number];
      const size_t grid_points_per_element = alg::accumulate(
          initial_extents[block.id()], 1_st, std::multiplies<size_t>());
      const auto initial_ref_levs = initial_refinement_levels[block.id()];
      const std::vector<ElementId<volume_dim>> element_ids =
          initial_element_ids_in_z_score_order(block.id(), initial_ref_levs);
      for (size_t i = 0; i < element_ids.size(); i++) {
        grid_points_by_element.push_back(grid_points_per_element);
        const auto& element_id = element_ids[i];
        const size_t target_proc =
            element_distribution.get_proc_for_element(element_id);
        dg_element_array(element_id)
            .insert(global_cache, initialization_items, target_proc);

        const size_t target_node =
            Parallel::node_of<size_t>(target_proc, local_cache);
        ++elements_per_core[target_proc];
        ++elements_per_node[target_node];
        grid_points_per_core[target_proc] += grid_points_per_element;
        grid_points_per_node[target_node] += grid_points_per_element;
        cost_per_core[target_proc] += cost_by_element_by_block[block_number][i];
        cost_per_node[target_node] += cost_by_element_by_block[block_number][i];
      }
    }
  } else {
    size_t which_proc = 0;
    for (const auto& block : domain.blocks()) {
      const size_t grid_points_per_element = alg::accumulate(
          initial_extents[block.id()], 1_st, std::multiplies<size_t>());
      const auto initial_ref_levs = initial_refinement_levels[block.id()];
      const std::vector<ElementId<volume_dim>> element_ids =
          initial_element_ids(block.id(), initial_ref_levs);
      while (procs_to_ignore.find(which_proc) != procs_to_ignore.end()) {
        which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
      }
      for (size_t i = 0; i < element_ids.size(); ++i) {
        dg_element_array(ElementId<volume_dim>(element_ids[i]))
            .insert(global_cache, initialization_items, which_proc);
        which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;

        const size_t target_node =
            Parallel::node_of<size_t>(which_proc, local_cache);
        ++elements_per_core[which_proc];
        ++elements_per_node[target_node];
        grid_points_per_core[which_proc] += grid_points_per_element;
        grid_points_per_node[target_node] += grid_points_per_element;
      }
    }
  }
  dg_element_array.doneInserting();

  if (use_z_order_distribution) {
    Parallel::printf(
        "\n%s\n", domain::diagnostic_info(
                      domain, local_cache, elements_per_core, elements_per_node,
                      grid_points_per_core, grid_points_per_node, cost_per_core,
                      cost_per_node));
  } else {
    Parallel::printf(
        "\n%s\n", domain::diagnostic_info(
                      domain, local_cache, elements_per_core, elements_per_node,
                      grid_points_per_core, grid_points_per_node));
  }
}
