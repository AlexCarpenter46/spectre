// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/Actions/LimitTimeStepLts.hpp"

#include <pup.h>

namespace control_system::StepChoosers {
PUP::able::PUP_ID LimitTimeStepLts::my_PUP_ID = 0;  // NOLINT
}  // namespace control_system::StepChoosers
