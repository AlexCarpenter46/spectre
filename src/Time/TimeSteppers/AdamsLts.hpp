// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/container/small_vector.hpp>
#include <cstddef>
#include <tuple>

#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsCoefficients.hpp"

/// \cond
namespace TimeSteppers {
template <typename T>
class BoundaryHistoryEvaluator;
class ConstBoundaryHistoryTimes;
class MutableBoundaryHistoryTimes;
}  // namespace TimeSteppers
namespace gsl {
template <class T>
class not_null;
}  // namespace gsl
/// \endcond

/// Shared LTS implementation for the two Adams-based methods.
namespace TimeSteppers::adams_lts {
/// Storage for LTS coefficients that should not allocate in typical
/// cases.
/// @{

// For order-k 2:1 stepping, in each small step, half the points will
// require interpolation (k entries each) and the others will not (1
// entry each).  So (2 small steps) * ((k/2 interpolations) * (from k
// points) + (k/2 non-interpolations)) = k (k + 1).  (The small steps
// round k/2 in different directions and the effect cancels out.)
constexpr size_t lts_coefficients_static_size =
    adams_coefficients::maximum_order * (adams_coefficients::maximum_order + 1);
using LtsCoefficients =
    boost::container::small_vector<std::tuple<TimeStepId, TimeStepId, double>,
                                   lts_coefficients_static_size>;
/// @}

/// Add the LTS boundary terms for to \p result for the given set of
/// coefficients.
template <typename T>
void apply_coefficients(gsl::not_null<T*> result,
                        const LtsCoefficients& coefficients,
                        const BoundaryHistoryEvaluator<T>& coupling);

/// Type of step to calculate coefficients for.
enum class StepType { Explicit, Predictor, Corrector };

#if defined(__clang__) and __clang__ < 17
#pragma GCC diagnostic push
// Clang is upset about \tparam for some reason.  Doxygen is happy.
#pragma GCC diagnostic ignored "-Wdocumentation"
#endif
/// Calculate the nonzero terms in the LTS boundary contribution.
///
/// All entries will have distinct pairs of IDs.
///
/// \tparam TimeType The type `Time` for a substep or
/// `ApproximateTime` for dense output.
///
/// \returns Container of tuples of `(local_id, remote_id, coefficient)`.
template <typename TimeType>
LtsCoefficients lts_coefficients(const ConstBoundaryHistoryTimes& local_times,
                                 const ConstBoundaryHistoryTimes& remote_times,
                                 const TimeType& end_time, StepType step_type);
#if defined(__clang__) and __clang__ < 17
#pragma GCC diagnostic pop
#endif

/// Remove old entries from the BoundaryHistory.
///
/// Removes any entries that are unnecessary for taking explicit steps
/// of order \p integration_order starting from the end of \p
/// local_times.  The order is not taken from the history because it
/// will be lower for predictor phases.
///
/// For a predictor-corrector method the corrector needs a superset of
/// the data used by the predictor, so no cleaning is needed in the
/// corrector phase.
void clean_boundary_history(const MutableBoundaryHistoryTimes& local_times,
                            const MutableBoundaryHistoryTimes& remote_times,
                            const size_t integration_order);
}  // namespace TimeSteppers::adams_lts
