// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

/// \cond
class LtsTimeStepper;
/// \endcond

namespace TimeStepperTestUtils::lts {
namespace patterns {
struct Lts2to1;
struct Lts3and1to2;
}  // namespace patterns

// FIXME doc all
template <typename Pattern>
void test_convergence(const LtsTimeStepper& stepper,
                      const std::pair<int32_t, int32_t>& step_range,
                      int32_t stride, bool output = false);

/// Test boundary computations with the same step size on both
/// neighbors.
void test_equal_rate(const LtsTimeStepper& stepper, size_t order,
                     size_t number_of_past_steps, double epsilon, bool forward);

/// Test the accuracy of boundary dense output.
void test_dense_output(const LtsTimeStepper& stepper);
}  // namespace TimeStepperTestUtils::lts
