// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "Evolution/Systems/Cce/NewmanPenrose.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"

namespace Cce {
namespace {
void pypp_test_volume_weyl() {
  pypp::SetupLocalPythonEnvironment local_python_env{"Evolution/Systems/Cce/"};

  const size_t num_pts = 5;

  pypp::check_with_random_values<1>(&(VolumeWeyl<Tags::Psi0>::apply),
                                    "NewmanPenrose", {"psi0"}, {{{1.0, 5.0}}},
                                    DataVector{num_pts});
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.Cce.NewmanPenrose", "[Unit][Cce]") {
  pypp_test_volume_weyl();
}

}  // namespace Cce
