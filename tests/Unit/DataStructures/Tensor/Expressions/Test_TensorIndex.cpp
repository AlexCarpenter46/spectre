// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/Tensor/Expressions/TensorIndex.hpp"
#include "DataStructures/Tensor/Expressions/TimeIndex.hpp"

SPECTRE_TEST_CASE("Unit.DataStructures.Tensor.Expression.TensorIndex",
                  "[DataStructures][Unit]") {
  // Test `make_tensorindex_list`
  // Check at compile time since some other tests use this metafunction
  static_assert(
      std::is_same_v<
          make_tensorindex_list<ti_j, ti_A, ti_b>,
          tmpl::list<std::decay_t<decltype(ti_j)>, std::decay_t<decltype(ti_A)>,
                     std::decay_t<decltype(ti_b)>>>,
      "make_tensorindex_list failed for non-empty list");
  static_assert(std::is_same_v<make_tensorindex_list<>, tmpl::list<>>,
                "make_tensorindex_list failed for empty list");

  // Test `get_tensorindex_value_with_opposite_valence`
  //
  // For (1) lower spacetime indices, (2) upper spacetime indices, (3) lower
  // spatial indices, and (4) upper spatial indices, the encoding of the index
  // with the same index type but opposite valence is checked for the following
  // cases: (i) the smallest encoding, (ii) the largest encoding, and (iii) a
  // value in between.

  // Lower spacetime
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(0) ==
        tenex::TensorIndex_detail::upper_sentinel);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::upper_sentinel - 1) ==
        tenex::TensorIndex_detail::spatial_sentinel - 1);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(115) ==
        tenex::TensorIndex_detail::upper_sentinel + 115);

  // Upper spacetime
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::upper_sentinel) == 0);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::spatial_sentinel - 1) ==
        tenex::TensorIndex_detail::upper_sentinel - 1);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::upper_sentinel + 88) == 88);

  // Lower spatial
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::spatial_sentinel) ==
        tenex::TensorIndex_detail::upper_spatial_sentinel);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::upper_spatial_sentinel - 1) ==
        tenex::TensorIndex_detail::max_sentinel - 1);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::spatial_sentinel + 232) ==
        tenex::TensorIndex_detail::upper_spatial_sentinel + 232);

  // Upper spatial
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::upper_spatial_sentinel) ==
        tenex::TensorIndex_detail::spatial_sentinel);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::max_sentinel - 1) ==
        tenex::TensorIndex_detail::upper_spatial_sentinel - 1);
  CHECK(tenex::get_tensorindex_value_with_opposite_valence(
            tenex::TensorIndex_detail::upper_spatial_sentinel + 3) ==
        tenex::TensorIndex_detail::spatial_sentinel + 3);

  // Test tensorindex_list_is_valid
  CHECK(tenex::tensorindex_list_is_valid<make_tensorindex_list<>>::value);
  CHECK(tenex::tensorindex_list_is_valid<make_tensorindex_list<ti_J>>::value);
  CHECK(tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_a, ti_c, ti_I, ti_B>>::value);
  CHECK(tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_t, ti_T, ti_T, ti_T, ti_t>>::value);
  CHECK(tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_d, ti_T, ti_D>>::value);
  CHECK(not tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_I, ti_a, ti_I>>::value);

  // Test tensorindex_list_is_valid
  CHECK(tenex::tensorindex_list_is_valid<make_tensorindex_list<>>::value);
  CHECK(tenex::tensorindex_list_is_valid<make_tensorindex_list<ti_J>>::value);
  CHECK(tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_a, ti_c, ti_I, ti_B>>::value);
  CHECK(tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_t, ti_T, ti_T, ti_T, ti_t>>::value);
  CHECK(tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_d, ti_T, ti_D>>::value);
  CHECK(not tenex::tensorindex_list_is_valid<
        make_tensorindex_list<ti_I, ti_a, ti_I>>::value);

  // Test generic_indices_at_same_positions
  CHECK(
      tenex::generic_indices_at_same_positions<make_tensorindex_list<>,
                                               make_tensorindex_list<>>::value);
  CHECK(tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_a, ti_c, ti_I, ti_B>,
        make_tensorindex_list<ti_a, ti_c, ti_I, ti_B>>::value);
  CHECK(not tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_a, ti_c, ti_I, ti_B>,
        make_tensorindex_list<ti_a, ti_c, ti_i, ti_B>>::value);
  CHECK(not tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_a, ti_c, ti_I, ti_B>,
        make_tensorindex_list<ti_a, ti_c, ti_I>>::value);
  CHECK(not tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_a, ti_c, ti_I>,
        make_tensorindex_list<ti_a, ti_c, ti_I, ti_B>>::value);
  CHECK(not tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_j, ti_B>,
        make_tensorindex_list<ti_B, ti_j>>::value);
  CHECK(tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_T>, make_tensorindex_list<ti_t>>::value);
  CHECK(tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_t, ti_T, ti_T>,
        make_tensorindex_list<ti_T, ti_t, ti_T>>::value);
  CHECK(tenex::generic_indices_at_same_positions<
        make_tensorindex_list<ti_i, ti_t, ti_C>,
        make_tensorindex_list<ti_i, ti_T, ti_C>>::value);
}
