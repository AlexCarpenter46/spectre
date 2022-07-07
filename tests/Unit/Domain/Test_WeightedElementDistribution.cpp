// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <set>
#include <unordered_set>
#include <vector>

#include "Domain/Block.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/Brick.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementDistribution.hpp"
#include "Domain/Structure/CreateInitialMesh.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Domain/WeightedElementDistribution.hpp"
#include "Domain/ZCurveIndex.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/Rational.hpp"

namespace {
std::vector<std::vector<double>> get_uniform_cost(const size_t num_blocks,
                                                  const size_t num_elements,
                                                  const double cost) {
  std::vector<std::vector<double>> costs_by_element_by_block(num_blocks);
  for (size_t i = 0; i < num_blocks; i++) {
    costs_by_element_by_block[i] = std::vector<double>(num_elements);
    for (size_t j = 0; j < num_elements; j++) {
      costs_by_element_by_block[i][j] = cost;
    }
  }
  return costs_by_element_by_block;
}

void print_costs_by_element_by_block(
    const std::vector<std::vector<double>>& costs_by_element_by_block) {
  std::cout << "Costs by element by block" << std::endl;
  for (size_t i = 0; i < costs_by_element_by_block.size(); i++) {
    std::cout << "Block " << i << ":\n\t{" << costs_by_element_by_block[i][0];
    for (size_t j = 1; j < costs_by_element_by_block[i].size(); j++) {
      std::cout << ", " << costs_by_element_by_block[i][j];
    }
    std::cout << "}" << std::endl;
  }
  std::cout << std::endl;
}

// std::vector<std::vector<double> > get_test_cost_2() {}

// std::vector<std::vector<double> > get_test_cost_3() {}

void print_element_distribution(
    const std::vector<std::vector<std::pair<size_t, size_t>>>&
        block_element_distribution) {
  std::cout << "Block element distribution" << std::endl;
  for (size_t i = 0; i < block_element_distribution.size(); i++) {
    std::cout << "Block " << i << ":\n\t{{"
              << block_element_distribution[i][0].first << ", "
              << block_element_distribution[i][0].second << "}";
    for (size_t j = 1; j < block_element_distribution[i].size(); j++) {
      std::cout << ", {" << block_element_distribution[i][j].first << ", "
                << block_element_distribution[i][j].second << "}";
    }
    std::cout << "}" << std::endl;
  }
  std::cout << std::endl;
}

template <size_t Dim>
void test(const size_t number_of_procs_with_elements,
          const std::vector<std::vector<double>>& cost_by_element_by_block,
          const std::unordered_set<size_t>& global_procs_to_ignore = {}) {
  // print_costs_by_element_by_block(cost_by_element_by_block);
  const domain::WeightedBlockZCurveProcDistribution<Dim> element_distribution{
      number_of_procs_with_elements, cost_by_element_by_block,
      global_procs_to_ignore};
  // // std::cout << element_distribution.block_element_distribution() <<
  // // std::endl;
  // print_element_distribution(element_distribution.block_element_distribution());
}

// for Sphere
template <size_t Dim>
std::vector<std::array<size_t, Dim>> get_initial_refinement_levels(
    const std::vector<Block<Dim>>& blocks, const size_t initial_refinement) {
  const size_t num_blocks = blocks.size();
  std::vector<std::array<size_t, Dim>> initial_refinement_levels(num_blocks);
  std::array<size_t, Dim> refinement_levels{};
  std::fill(refinement_levels.begin(), refinement_levels.end(),
            initial_refinement);

  for (size_t i = 0; i < num_blocks; i++) {
    initial_refinement_levels[i] = refinement_levels;
  }

  return initial_refinement_levels;
}

// for Brick
template <size_t Dim>
std::vector<std::array<size_t, Dim>> get_initial_refinement_levels(
    const std::vector<Block<Dim>>& blocks,
    const std::array<size_t, 3>& initial_refinement_level_xyz) {
  const size_t num_blocks = blocks.size();
  std::vector<std::array<size_t, Dim>> initial_refinement_levels(num_blocks);

  for (size_t i = 0; i < num_blocks; i++) {
    initial_refinement_levels[i] = initial_refinement_level_xyz;
  }

  return initial_refinement_levels;
}

template <size_t Dim, bool IsWeighted = true,
          typename ElementDistribution = tmpl::conditional_t<
              IsWeighted, domain::WeightedBlockZCurveProcDistribution<Dim>,
              domain::BlockZCurveProcDistribution<Dim>>>
ElementDistribution get_element_distribution(
    const std::vector<std::array<size_t, Dim>>& initial_refinement_levels,
    const size_t num_of_procs_to_use,
    const std::unordered_set<size_t>& procs_to_ignore) {
  return ElementDistribution{num_of_procs_to_use, initial_refinement_levels,
                             procs_to_ignore};
}

template <size_t Dim, bool IsWeighted = true>
void test_z_curve_index(
    const Domain<Dim>& domain,
    const std::array<size_t, 3>& initial_refinement_level_xyz,
    const size_t num_of_procs_to_use,
    const std::unordered_set<size_t>& procs_to_ignore) {
  const std::vector<Block<Dim>>& blocks = domain.blocks();

  const std::vector<std::array<size_t, Dim>> initial_refinement_levels =
      get_initial_refinement_levels(blocks, initial_refinement_level_xyz);

  // size_t num_elements = two_to_the(initial_refinement_levels[0]);
  // for (size_t i = 1; i < Dim; i++) {
  //   num_elements *= two_to_the(initial_refinement_levels[1]);
  // }

  size_t num_elements = 0;
  for (size_t i = 0; i < initial_refinement_levels.size(); i++) {
    size_t num_elements_this_block =
        two_to_the(initial_refinement_levels[i][0]);
    for (size_t j = 1; j < Dim; j++) {
      num_elements_this_block *= two_to_the(initial_refinement_levels[i][j]);
    }
    num_elements += num_elements_this_block;
  }

  const auto element_distribution = get_element_distribution<Dim, IsWeighted>(
      initial_refinement_levels, num_of_procs_to_use, procs_to_ignore);

  // size_t run = 0;
  for (const auto& block : blocks) {
    // std::cout << "block id : " << block.id() << std::endl;

    const auto initial_ref_levs = initial_refinement_levels[block.id()];
    const std::vector<ElementId<Dim>> element_ids =
        initial_element_ids(block.id(), initial_ref_levs);
    // for (const auto& element_id : element_ids) {
    for (size_t j = 0; j < num_elements; j++) {
      const auto& element_id = element_ids[j];
      // std::cout << "element_id : " << element_id << std::endl;

      // const size_t target_proc =
      //     element_distribution.get_proc_for_element(element_id);

      const size_t result_z_order_index =
          domain::z_curve_index_from_element_id(element_id);
      // std::cout << "result_z_order_index : " << result_z_order_index <<
      // std::endl;
      const std::array<size_t, Dim> result_element_id =
          domain::element_id_from_z_curve_index(result_z_order_index,
                                                initial_ref_levs);

      std::array<size_t, Dim> expected_element_id;
      for (size_t i = 0; i < Dim; ++i) {
        expected_element_id[i] = element_id.segment_id(i).index();
      }

      CHECK(result_element_id == expected_element_id);

      // // std::cout << "target_proc : " << target_proc << std::endl;
      // if (run > 1) break;
      // run++;
      // std::cout << std::endl;
    }
    const std::vector<ElementId<Dim>> element_ids_in_z_score_order =
        initial_element_ids_in_z_score_order(block.id(), initial_ref_levs);
    for (size_t j = 0; j < num_elements; j++) {
      const auto& element_id = element_ids_in_z_score_order[j];
      CHECK(domain::z_curve_index_from_element_id(element_id) == j);
    }
    // std::cout << std::endl;
  }
}

template <size_t Dim>
void test_compute_minimum_grid_spacing(
    const std::vector<std::array<size_t, Dim>>& initial_extents,
    const std::vector<std::array<size_t, Dim>>& initial_refinement,
    const Domain<Dim>& domain, const Spectral::Quadrature& quadrature,
    const size_t array_index) {
  // ===== From DgDomain.hpp =====
  const ElementId<Dim> element_id{array_index};
  const auto& my_block = domain.blocks()[element_id.block_id()];
  Mesh<Dim> mesh = ::domain::Initialization::create_initial_mesh(
      initial_extents, element_id, quadrature);
  Element<Dim> element = ::domain::Initialization::create_initial_element(
      element_id, my_block, initial_refinement);
  ElementMap<Dim, Frame::Grid> element_map{
      element_id, my_block.is_time_dependent()
                      ? my_block.moving_mesh_logical_to_grid_map().get_clone()
                      : my_block.stationary_map().get_to_grid_frame()};

  std::unique_ptr<
      ::domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>>
      grid_to_inertial_map;
  if (my_block.is_time_dependent()) {
    grid_to_inertial_map =
        my_block.moving_mesh_grid_to_inertial_map().get_clone();
  } else {
    grid_to_inertial_map =
        ::domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
            ::domain::CoordinateMaps::Identity<Dim>{});
  }
  // ===========================

  // Get logical coordinates, i.e. domain::tags::LogicalCoordinates<Dim>
  tnsr::I<DataVector, Dim, Frame::ElementLogical> logical_coords{};
  domain::Tags::LogicalCoordinates<Dim>::function(
      make_not_null(&logical_coords), mesh);

  // Get grid coordinates, i.e.
  //     domain::tags::MappedCoordinates<
  //         domain::Tags::ElementMap<Dim, Frame::Grid>,
  //         domain::Tags::Coordinates<Dim, Frame::ElementLogical>>
  //
  // (aka domain::Tags::Coordinates<Dim, Frame::Grid>)
  tnsr::I<DataVector, Dim, Frame::Grid> grid_coords{};
  domain::Tags::MappedCoordinates<
      domain::Tags::ElementMap<Dim, Frame::Grid>,
      domain::Tags::Coordinates<Dim, Frame::ElementLogical>>::
      function(make_not_null(&grid_coords), element_map, logical_coords);

  // Get ::Tags::Time

  // Get domain::Tags::FunctionsOfTime

  // Get domain::Tags::CoordinatesMeshVelocityAndJacobians<Dim>

  // Get inertial coordinates, i.e.
  //     domain::Tags::InertialFromGridCoordinatesCompute<Dim>
  //
  // (aka domain::Tags::Coordinates<Dim, Frame::Inertial>)

  // Get minimum grid spacing, i.e.
  //     domain::Tags::MinimumGridSpacingCompute<Dim, Frame::Inertial>
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.WeightedElementDistribution", "[Domain][Unit]") {
  // test<3>(10, get_uniform_cost(6, 4, 1.0));

  const size_t Dim = 3;

  // const double inner_radius = 10.0;
  // const double outer_radius = 110.0;
  // const size_t initial_refinement = 1;
  // const std::array<size_t, 2> initial_number_of_grid_points{{5, 7}};
  // const bool use_equiangular_map = false;

  // domain::creators::Sphere sphere(
  //     inner_radius, outer_radius, initial_refinement,
  //     initial_number_of_grid_points, use_equiangular_map);

  const std::array<double, 3> lower_xyz = {0.0, 0.0, 0.0};
  const std::array<double, 3> upper_xyz = {1.0, 10.0, 100.0};
  const std::array<size_t, 3> initial_refinement_level_xyz = {1, 2, 3};
  const std::array<size_t, 3> initial_number_of_grid_points_in_xyz = {2, 4, 6};
  const std::array<bool, 3> is_periodic_in_xyz = {{false, false, false}};

  domain::creators::Brick brick(
      lower_xyz, upper_xyz, initial_refinement_level_xyz,
      initial_number_of_grid_points_in_xyz, is_periodic_in_xyz);

  const size_t num_of_procs_to_use = 3;
  const std::unordered_set<size_t> procs_to_ignore{};
  // std::vector<std::array<size_t, Dim>> initial_refinement_levels

  // const domain::BlockZCurveProcDistribution<Dim> element_distribution{
  //     num_of_procs_to_use, initial_refinement_levels, procs_to_ignore};

  test_z_curve_index<Dim, false>(brick.create_domain(),
                                 initial_refinement_level_xyz,
                                 num_of_procs_to_use, procs_to_ignore);
}
