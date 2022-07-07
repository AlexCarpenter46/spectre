// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/ZCurveIndex.hpp"

#include <array>
#include <cstddef>
#include <utility>

#include "Domain/Structure/ElementId.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace domain {

template <size_t Dim>
size_t z_curve_index_from_element_id(const ElementId<Dim>& element_id) {
  // for the bit manipulation of the element index, we need to interleave the
  // indices in each dimension in order according to how many bits are in the
  // index representation. This variable stores the refinement level and
  // dimension index in ascending order of refinement level, representing a
  // permutation of the dimensions
  // pair<refinement level, dim index> in order of ascending refinement
  std::array<std::pair<size_t, size_t>, Dim>
      dimension_by_highest_refinement_level;
  for (size_t i = 0; i < Dim; ++i) {
    dimension_by_highest_refinement_level.at(i) =
        std::make_pair(element_id.segment_id(i).refinement_level(), i);
  }
  alg::sort(dimension_by_highest_refinement_level,
            [](const std::pair<size_t, size_t>& lhs,
               const std::pair<size_t, size_t>& rhs) {
              return lhs.first < rhs.first;
            });

  size_t element_order_index = 0;

  // 'gap' the lowest refinement direction bits as:
  // ... x1 x0 -> ... x1 0 0 x0,
  // then bitwise or in 'gap'ed and shifted next-lowest refinement direction
  // bits as:
  // ... y2 y1 y0 -> ... y2 0 y1 x1 0 y0 x0
  // then bitwise or in 'gap'ed and shifted highest-refinement direction bits
  // as:
  // ... z3 z2 z1 z0 -> z3 z2 y2 z1 y1 x1 z0 y0 x0
  // note that we must skip refinement-level 0 dimensions as though they are
  // not present
  size_t leading_gap = 0;
  for (size_t i = 0; i < Dim; ++i) {
    const size_t id_to_gap_and_shift =
        element_id
            .segment_id(
                gsl::at(dimension_by_highest_refinement_level, i).second)
            .index();
    size_t total_gap = leading_gap;
    if (gsl::at(dimension_by_highest_refinement_level, i).first > 0) {
      ++leading_gap;
    }
    for (size_t bit_index = 0;
         bit_index < gsl::at(dimension_by_highest_refinement_level, i).first;
         ++bit_index) {
      // This operation will not overflow for our present use of `ElementId`s.
      // This technique densely assigns an ElementID a unique size_t identifier
      // determining the Morton curve order, and `ElementId` supports refinement
      // levels such that a global index within a block will fit in a 64-bit
      // unsigned integer.
      element_order_index |=
          ((id_to_gap_and_shift & two_to_the(bit_index)) << total_gap);
      for (size_t j = 0; j < Dim; ++j) {
        if (i != j and
            bit_index + 1 <
                gsl::at(dimension_by_highest_refinement_level, j).first) {
          ++total_gap;
        }
      }
    }
  }
  return element_order_index;
}

template <size_t Dim>
std::array<size_t, Dim> element_id_from_z_curve_index(
    const size_t z_order_index,
    const std::array<size_t, Dim>& block_refinements) {
  std::array<std::pair<size_t, size_t>, Dim>
      dimension_by_highest_refinement_level;
  for (size_t i = 0; i < Dim; ++i) {
    dimension_by_highest_refinement_level.at(i) =
        std::make_pair(gsl::at(block_refinements, i), i);
  }
  // {{L, Dim}, {L, Dim}, ...}
  alg::sort(dimension_by_highest_refinement_level,
            [](const std::pair<size_t, size_t>& lhs,
               const std::pair<size_t, size_t>& rhs) {
              return lhs.first < rhs.first;
            });

  std::array<std::pair<size_t, size_t>, Dim>
      segment_indices_by_highest_refinement_level;
  for (size_t i = 0; i < Dim; ++i) {
    segment_indices_by_highest_refinement_level.at(i) = std::make_pair(
        0, gsl::at(dimension_by_highest_refinement_level, i).second);
  }
  size_t starting_dim_index = 0;
  size_t bit_index = 0;
  size_t element_order_index = z_order_index;

  // TODO : try to optimize by handling last case separately when it's the
  // remaining bits of the highest refinement, i.e. probably do:
  // starting_dim_index < (Dim - 1) then handle the last case separately after
  while (starting_dim_index < Dim) {
    const size_t refinement_level =
        gsl::at(dimension_by_highest_refinement_level, starting_dim_index)
            .first;
    // TODO : this is probably wrong
    // const size_t num_refinement_bits =
    //     bit_index == 0 ? two_to_the(refinement_level) :
    //         two_to_the(refinement_level) - two_to_the(bit_index);
    const size_t num_refinement_bits = refinement_level - bit_index;
    for (size_t i = 0; i < num_refinement_bits; i++) {
      for (size_t dim_index = starting_dim_index; dim_index < Dim;
           dim_index++) {
        segment_indices_by_highest_refinement_level.at(dim_index).first |=
            ((element_order_index & 1) << bit_index);
        element_order_index >>= 1;
      }
      bit_index++;
    }

    starting_dim_index++;
    while (starting_dim_index < Dim and
           refinement_level == gsl::at(dimension_by_highest_refinement_level,
                                       starting_dim_index)
                                   .first) {
      starting_dim_index++;
    }
  }

  std::array<size_t, Dim> element_id{};

  for (size_t i = 0; i < Dim; i++) {
    gsl::at(element_id,
            gsl::at(segment_indices_by_highest_refinement_level, i).second) =
        gsl::at(segment_indices_by_highest_refinement_level, i).first;
  }

  return element_id;
}
#define GET_DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                              \
  template size_t z_curve_index_from_element_id(                            \
      const ElementId<GET_DIM(data)>& element_id);                          \
  template std::array<size_t, GET_DIM(data)> element_id_from_z_curve_index( \
      const size_t z_order_index,                                           \
      const std::array<size_t, GET_DIM(data)>& block_refinements);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

}  // namespace domain
