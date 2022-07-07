// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/WeightedElementDistribution.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

#include "Domain/Structure/ElementId.hpp"
#include "Domain/ZCurveIndex.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace domain {
template <size_t Dim>
WeightedBlockZCurveProcDistribution<Dim>::WeightedBlockZCurveProcDistribution(
    const size_t number_of_procs_with_elements,
    const std::vector<std::vector<double> >& cost_by_element_by_block,
    const std::unordered_set<size_t>& global_procs_to_ignore) {
  block_element_distribution_ =
      std::vector<std::vector<std::pair<size_t, size_t> > >(
          cost_by_element_by_block.size());

  double total_cost = 0.0;
  double min_cost = std::numeric_limits<double>::max();
  double max_cost = std::numeric_limits<double>::min();

  for (auto& block : cost_by_element_by_block) {
    const auto min_element_this_block = *alg::min_element(block);
    const auto max_element_this_block = *alg::max_element(block);
    if (min_element_this_block < min_cost) {
      min_cost = min_element_this_block;
    }
    if (max_element_this_block > max_cost) {
      max_cost = max_element_this_block;
    }

    for (double element_cost : block) {
      total_cost += element_cost;
    }
  }

  size_t current_block = 0;
  size_t current_element_of_current_block = 0;
  double cost_remaining = total_cost;
  // This variable will keep track of how many global procs we've skipped over
  // so far. This bookkeeping is necessary so the element gets placed on the
  // correct global proc. The loop variable `i` does not correspond to global
  // proc number. It's just an index
  size_t number_of_ignored_procs_so_far = 0;
  for (size_t i = 0; i < number_of_procs_with_elements; ++i) {
    size_t global_proc_number = i + number_of_ignored_procs_so_far;
    while (global_procs_to_ignore.find(global_proc_number) !=
           global_procs_to_ignore.end()) {
      ++number_of_ignored_procs_so_far;
      ++global_proc_number;
    }

    // initialize cost for this proc to be the current element
    // double cost_spent_on_proc =
    //     cost_by_element_by_block[current_block]
    //         [current_element_of_current_block];
    // size_t num_elements_distributed_to_proc = 1;
    double target_cost_per_proc =
        cost_remaining / (number_of_procs_with_elements - i);
    double cost_spent_on_proc = 0.0;
    size_t total_elements_distributed_to_proc = 0;
    bool add_more_elements_to_proc = true;
    const size_t num_blocks = cost_by_element_by_block.size();
    // while we still have cost allowed on the proc
    while (add_more_elements_to_proc and (current_block < num_blocks)) {
      const size_t num_elements_current_block =
          cost_by_element_by_block[current_block].size();
      // while we still have elements left on the block and we still
      // have cost allowed on the proc
      size_t num_elements_distributed_to_proc = 0;

      while (add_more_elements_to_proc and
             (current_element_of_current_block < num_elements_current_block)) {
        const double element_cost =
            cost_by_element_by_block[current_block]
                                    [current_element_of_current_block];

        if (total_elements_distributed_to_proc == 0) {
          cost_remaining -= element_cost;
          cost_spent_on_proc = element_cost;
          num_elements_distributed_to_proc = 1;
          total_elements_distributed_to_proc = 1;
          current_element_of_current_block++;
        } else {
          const double current_cost_diff =
              abs(target_cost_per_proc - cost_spent_on_proc);
          const double next_cost_diff =
              abs(target_cost_per_proc - (cost_spent_on_proc + element_cost));

          if (current_cost_diff <= next_cost_diff) {
            add_more_elements_to_proc = false;
          } else {
            cost_spent_on_proc += element_cost;
            cost_remaining -= element_cost;
            num_elements_distributed_to_proc++;
            total_elements_distributed_to_proc++;
            current_element_of_current_block++;
          }
        }
      }

      block_element_distribution_.at(current_block)
          .emplace_back(std::make_pair(global_proc_number,
                                       num_elements_distributed_to_proc));
      if (current_element_of_current_block >= num_elements_current_block) {
        ++current_block;
        current_element_of_current_block = 0;
      }
    }
  }
}

template <size_t Dim>
size_t WeightedBlockZCurveProcDistribution<Dim>::get_proc_for_element(
    const ElementId<Dim>& element_id) const {
  const size_t element_order_index = z_curve_index_from_element_id(element_id);
  size_t total_so_far = 0;
  for (const std::pair<size_t, size_t>& element_info :
       gsl::at(block_element_distribution_, element_id.block_id())) {
    if (total_so_far <= element_order_index and
        element_info.second + total_so_far > element_order_index) {
      return element_info.first;
    }
    total_so_far += element_info.second;
  }
  ERROR(
      "Processor not successfully chosen. This indicates a flaw in the logic "
      "of WeightedBlockZCurveProcDistribution.");
}

#define GET_DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) \
  template class WeightedBlockZCurveProcDistribution<GET_DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef GET_DIM
#undef INSTANTIATION
}  // namespace domain
