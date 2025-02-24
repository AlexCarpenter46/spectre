// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "DataStructures/Tensor/TypeAliases.hpp"

/// \cond
namespace ylm {
template <typename Frame>
class Strahlkorper;
}  // namespace ylm
template <size_t Dim>
class Domain;
namespace domain::FunctionsOfTime {
class FunctionOfTime;
}  // namespace domain::FunctionsOfTime
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
/// \endcond

// Todo: Add brief and better comment.
// Transforms cartesian coordinates of a strahlkorper
// from one frame to another, by calling block_logical_coordinates and
// calling the correct map functions.
template <typename SrcFrame, typename DestFrame>
void coords_to_different_frame(
    gsl::not_null<tnsr::I<DataVector, 3, DestFrame>*> dest_cartesian_coords,
    const tnsr::I<DataVector, 3, SrcFrame>& src_cartesian_coords,
    const Domain<3>& domain,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time,
    double time);
