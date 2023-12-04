// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>

/// \cond
namespace amr {
template <size_t VolumeDim>
struct Info;
}  // namespace amr
template <size_t VolumeDim>
class Direction;
template <size_t VolumeDim>
class ElementId;
template <size_t VolumeDim>
class Mesh;
template <size_t VolumeDim>
class Neighbors;
/// \endcond

namespace amr {
/// \ingroup AmrGroup
/// \brief returns the ElementId and Mesh of the new neighbors in the given
/// `direction` of the Element whose ElementId is `my_id` given the
/// `previous_neighbors_in_direction` and their amr::Info\.
///
/// \note `previous_neighbors_in_direction` should be from the parent (or a
/// child) of the Element with `my_id` if `my_id` corresponds to a newly created
/// child (or parent) Element.
///
/// \note `previous_neighbors_amr_info` may contain info from neighbors in
/// directions other than `direction`
template <size_t VolumeDim>
std::unordered_map<ElementId<VolumeDim>, Mesh<VolumeDim>> new_neighbor_ids(
    const ElementId<VolumeDim>& my_id, const Direction<VolumeDim>& direction,
    const Neighbors<VolumeDim>& previous_neighbors_in_direction,
    const std::unordered_map<ElementId<VolumeDim>, Info<VolumeDim>>&
        previous_neighbors_amr_info);
}  // namespace amr
