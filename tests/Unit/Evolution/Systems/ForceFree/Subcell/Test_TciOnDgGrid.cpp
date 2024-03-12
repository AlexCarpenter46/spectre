// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <limits>
#include <memory>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DgSubcell/Mesh.hpp"
#include "Evolution/DgSubcell/Projection.hpp"
#include "Evolution/DgSubcell/SubcellOptions.hpp"
#include "Evolution/DgSubcell/Tags/Mesh.hpp"
#include "Evolution/DgSubcell/Tags/SubcellOptions.hpp"
#include "Evolution/Systems/ForceFree/Subcell/TciOnDgGrid.hpp"
#include "Evolution/Systems/ForceFree/Subcell/TciOptions.hpp"
#include "Evolution/Systems/ForceFree/System.hpp"
#include "Evolution/Systems/ForceFree/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace {
enum class TestThis {
  AllGood,
  PerssonMagTildeE,
  PerssonMagTildeB,
  PerssonTildeQ,
  PerssonTildeQBelowCutoff,
  PerssonTildeQDoNotCheck,
  RdmpMagTildeE,
  RdmpMagTildeB
};

void test(const TestThis test_this, const int expected_tci_status) {
  CAPTURE(test_this);
  CAPTURE(expected_tci_status);

  const Mesh<3> dg_mesh{6, Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};
  const Mesh<3> subcell_mesh = evolution::dg::subcell::fd::mesh(dg_mesh);

  using TildeE = ForceFree::Tags::TildeE;
  using TildeB = ForceFree::Tags::TildeB;
  using TildeQ = ForceFree::Tags::TildeQ;

  using VarsForTciTest = Variables<tmpl::list<TildeE, TildeB, TildeQ>>;
  VarsForTciTest dg_vars{dg_mesh.number_of_grid_points(), 0.0};

  // set variables on the dg mesh for the test
  get(get<TildeQ>(dg_vars)) = 1.0;
  for (size_t i = 0; i < 3; ++i) {
    get<TildeE>(dg_vars).get(i) = 2.0;
    get<TildeB>(dg_vars).get(i) = 3.0;
  }

  // Just use flat spacetime since none of the TCI checks depends on metric
  // variables
  tnsr::ii<DataVector, 3, Frame::Inertial> spatial_metric{
      dg_mesh.number_of_grid_points(), 0.0};
  tnsr::II<DataVector, 3, Frame::Inertial> inv_spatial_metric{
      dg_mesh.number_of_grid_points(), 0.0};
  for (size_t i = 0; i < 3; ++i) {
    spatial_metric.get(i, i) = inv_spatial_metric.get(i, i) = 1.0;
  }
  const Scalar<DataVector> sqrt_det_spatial_metric{
      dg_mesh.number_of_grid_points(), 1.0};

  const double persson_exponent = 4.0;

  const evolution::dg::subcell::SubcellOptions subcell_options{
      persson_exponent,
      1_st,
      1.0e-10,
      1.0e-10,
      false,
      evolution::dg::subcell::fd::ReconstructionMethod::DimByDim,
      false,
      std::nullopt,
      fd::DerivativeOrder::Two,
      1,
      1,
      1};

  const ForceFree::subcell::TciOptions tci_options{1.0e-10};

  auto box = db::create<db::AddSimpleTags<
      ::Tags::Variables<VarsForTciTest::tags_list>, ::domain::Tags::Mesh<3>,
      ::evolution::dg::subcell::Tags::Mesh<3>,
      gr::Tags::SqrtDetSpatialMetric<DataVector>,
      gr::Tags::SpatialMetric<DataVector, 3>,
      gr::Tags::InverseSpatialMetric<DataVector, 3>,
      ForceFree::subcell::Tags::TciOptions,
      evolution::dg::subcell::Tags::SubcellOptions<3>,
      evolution::dg::subcell::Tags::DataForRdmpTci>>(
      dg_vars, dg_mesh, subcell_mesh, sqrt_det_spatial_metric, spatial_metric,
      inv_spatial_metric, tci_options, subcell_options,
      evolution::dg::subcell::RdmpTciData{});

  const size_t point_to_change = dg_mesh.number_of_grid_points() / 2;

  if (test_this == TestThis::PerssonMagTildeE) {
    db::mutate<TildeE>(
        [point_to_change](const auto tilde_e_ptr) {
          (*tilde_e_ptr).get(0)[point_to_change] *= 2.0;
        },
        make_not_null(&box));
  } else if (test_this == TestThis::PerssonMagTildeB) {
    db::mutate<TildeB>(
        [point_to_change](const auto tilde_b_ptr) {
          (*tilde_b_ptr).get(0)[point_to_change] *= 2.0;
        },
        make_not_null(&box));
  } else if (test_this == TestThis::PerssonTildeQ) {
    db::mutate<TildeQ>(
        [point_to_change](const auto tilde_q_ptr) {
          get(*tilde_q_ptr)[point_to_change] *= 2.0;
        },
        make_not_null(&box));
  } else if (test_this == TestThis::PerssonTildeQBelowCutoff) {
    db::mutate<TildeQ>(
        [point_to_change](const auto tilde_q_ptr) {
          get(*tilde_q_ptr)[point_to_change] *= 2.0;
          get(*tilde_q_ptr) *= 1e-20;
        },
        make_not_null(&box));
  } else if (test_this == TestThis::PerssonTildeQDoNotCheck) {
    db::mutate<TildeQ, ForceFree::subcell::Tags::TciOptions>(
        [point_to_change](const auto tilde_q_ptr, const auto tci_options_ptr) {
          get(*tilde_q_ptr)[point_to_change] *= 2.0;
          tci_options_ptr->tilde_q_cutoff = std::nullopt;
        },
        make_not_null(&box));
  }

  // Set the RDMP TCI past data.
  using std::max;
  using std::min;

  const auto dg_mag_tilde_e = get(magnitude(db::get<TildeE>(box)));
  const auto dg_mag_tilde_b = get(magnitude(db::get<TildeB>(box)));
  const auto dg_tilde_q = get(db::get<TildeQ>(box));
  const auto subcell_mag_tilde_e = evolution::dg::subcell::fd::project(
      dg_mag_tilde_e, dg_mesh, subcell_mesh.extents());
  const auto subcell_mag_tilde_b = evolution::dg::subcell::fd::project(
      dg_mag_tilde_b, dg_mesh, subcell_mesh.extents());

  evolution::dg::subcell::RdmpTciData past_rdmp_tci_data{};
  past_rdmp_tci_data.max_variables_values =
      DataVector{max(max(dg_mag_tilde_e), max(subcell_mag_tilde_e)),
                 max(max(dg_mag_tilde_b), max(subcell_mag_tilde_b))};
  past_rdmp_tci_data.min_variables_values =
      DataVector{min(min(dg_mag_tilde_e), min(subcell_mag_tilde_e)),
                 min(min(dg_mag_tilde_b), min(subcell_mag_tilde_b))};

  const evolution::dg::subcell::RdmpTciData expected_rdmp_tci_data =
      past_rdmp_tci_data;

  // Modify past data if we are expecting an RDMP TCI failure.
  db::mutate<evolution::dg::subcell::Tags::DataForRdmpTci>(
      [&past_rdmp_tci_data, &test_this](const auto rdmp_tci_data_ptr) {
        *rdmp_tci_data_ptr = past_rdmp_tci_data;
        // Assumes min is positive, increase it so we fail the TCI
        if (test_this == TestThis::RdmpMagTildeE) {
          rdmp_tci_data_ptr->min_variables_values[0] *= 1.01;
        } else if (test_this == TestThis::RdmpMagTildeB) {
          rdmp_tci_data_ptr->min_variables_values[1] *= 1.01;
        }
      },
      make_not_null(&box));

  const bool element_stays_on_dg = false;
  const std::tuple<int, evolution::dg::subcell::RdmpTciData> result =
      db::mutate_apply<ForceFree::subcell::TciOnDgGrid>(
          make_not_null(&box), persson_exponent, element_stays_on_dg);

  CHECK(get<1>(result) == expected_rdmp_tci_data);

  if (test_this == TestThis::AllGood) {
    CHECK_FALSE(get<0>(result));
  } else {
    CHECK(get<0>(result) == expected_tci_status);
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.ForceFree.Systems.Subcell.TciOnDgGrid",
                  "[Unit][Evolution]") {
  test(TestThis::AllGood, 0);
  test(TestThis::PerssonMagTildeE, -1);
  test(TestThis::PerssonMagTildeB, -2);
  test(TestThis::PerssonTildeQ, -3);
  test(TestThis::PerssonTildeQBelowCutoff, 0);
  test(TestThis::PerssonTildeQDoNotCheck, 0);
  test(TestThis::RdmpMagTildeE, -4);
  test(TestThis::RdmpMagTildeB, -5);
}
