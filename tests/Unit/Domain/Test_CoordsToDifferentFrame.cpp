// // Distributed under the MIT License.
// // See LICENSE.txt for details.

// #include "Framework/TestingFramework.hpp"

// #include <array>
// #include <random>

// #include "DataStructures/DataVector.hpp"
// #include "DataStructures/Tensor/IndexType.hpp"
// #include "DataStructures/Tensor/Tensor.hpp"
// #include "DataStructures/Tensor/TypeAliases.hpp"
// #include "Domain/Creators/RegisterDerivedWithCharm.hpp"
// #include "Domain/Creators/Sphere.hpp"
// #include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
// #include "Domain/Creators/TimeDependence/Shape.hpp"
// #include "Domain/Domain.hpp"
// #include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
// #include "Domain/StrahlkorperTransformations.hpp"
// #include "Framework/TestHelpers.hpp"
// #include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
// #include "NumericalAlgorithms/SphericalHarmonics/StrahlkorperFunctions.hpp"
// #include
// "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrHorizon.hpp"
// #include "Utilities/Gsl.hpp"

// namespace {
// template <typename SrcFrame, typename DestFrame>
// void test_coords_to_different_frame() {

// }
// SPECTRE_TEST_CASE("Unit.Domain.StrahlkorperTransformations", "[Unit]") {
//   domain::creators::register_derived_with_charm();
//   domain::creators::time_dependence::register_derived_with_charm();
//   domain::FunctionsOfTime::register_derived_with_charm();
//   test_coords_to_different_frame<Frame::Grid, Frame::Inertial>();
//   test_coords_to_different_frame<Frame::Inertial, Frame::Distorted>();
//   test_coords_to_different_frame<Frame::Inertial, Frame::Grid>();
//   test_coords_to_different_frame<Frame::Grid, Frame::Inertial>();
//   test_coords_to_different_frame<Frame::Grid, Frame::Distorted>();
// }
// } // namespace
