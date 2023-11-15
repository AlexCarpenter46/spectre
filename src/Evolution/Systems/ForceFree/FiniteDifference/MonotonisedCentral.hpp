// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <pup.h>
#include <utility>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/Structure/DirectionIdMap.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DgSubcell/GhostData.hpp"
#include "Evolution/DgSubcell/Tags/GhostDataForReconstruction.hpp"
#include "Evolution/DgSubcell/Tags/Mesh.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/NormalCovectorAndMagnitude.hpp"
#include "Evolution/Systems/ForceFree/FiniteDifference/Reconstructor.hpp"
#include "Evolution/Systems/ForceFree/FiniteDifference/Tags.hpp"
#include "Evolution/Systems/ForceFree/Tags.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
template <size_t dim>
class Direction;
template <size_t dim>
class Element;
template <size_t dim>
class ElementId;
template <size_t dim>
class Mesh;
template <typename recons_tags>
class Variables;
namespace gsl {
template <typename>
class not_null;
}  // namespace gsl
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace ForceFree::fd {

/*!
 * \brief Monotonised central reconstruction. See
 * ::fd::reconstruction::monotonised_central() for details.
 */

class MonotonisedCentral : public Reconstructor {
 private:
  using TildeE = ForceFree::Tags::TildeE;
  using TildeB = ForceFree::Tags::TildeB;
  using TildePsi = ForceFree::Tags::TildePsi;
  using TildePhi = ForceFree::Tags::TildePhi;
  using TildeQ = ForceFree::Tags::TildeQ;
  using TildeJ = ForceFree::Tags::TildeJ;

  using volume_vars_tags =
      tmpl::list<TildeE, TildeB, TildePsi, TildePhi, TildeQ>;

  using recons_tags = ForceFree::fd::tags_list_for_reconstruction;

 public:
  static constexpr size_t dim = 3;

  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Monotonised central reconstruction scheme."};

  MonotonisedCentral() = default;
  MonotonisedCentral(MonotonisedCentral&&) = default;
  MonotonisedCentral& operator=(MonotonisedCentral&&) = default;
  MonotonisedCentral(const MonotonisedCentral&) = default;
  MonotonisedCentral& operator=(const MonotonisedCentral&) = default;
  ~MonotonisedCentral() override = default;

  explicit MonotonisedCentral(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(Reconstructor, MonotonisedCentral);

  auto get_clone() const -> std::unique_ptr<Reconstructor> override;

  static constexpr bool use_adaptive_order = false;

  void pup(PUP::er& p) override;

  size_t ghost_zone_size() const override { return 2; }

  using reconstruction_argument_tags =
      tmpl::list<::Tags::Variables<volume_vars_tags>, TildeJ,
                 domain::Tags::Element<dim>,
                 evolution::dg::subcell::Tags::GhostDataForReconstruction<dim>,
                 evolution::dg::subcell::Tags::Mesh<dim>>;

  void reconstruct(
      gsl::not_null<std::array<Variables<recons_tags>, dim>*>
          vars_on_lower_face,
      gsl::not_null<std::array<Variables<recons_tags>, dim>*>
          vars_on_upper_face,
      const Variables<volume_vars_tags>& volume_vars,
      const tnsr::I<DataVector, 3, Frame::Inertial>& tilde_j,
      const Element<dim>& element,
      const DirectionIdMap<dim, evolution::dg::subcell::GhostData>& ghost_data,
      const Mesh<dim>& subcell_mesh) const;

  void reconstruct_fd_neighbor(
      gsl::not_null<Variables<recons_tags>*> vars_on_face,
      const Variables<volume_vars_tags>& volume_vars,
      const tnsr::I<DataVector, 3, Frame::Inertial>& tilde_j,
      const Element<dim>& element,
      const DirectionIdMap<dim, evolution::dg::subcell::GhostData>& ghost_data,
      const Mesh<dim>& subcell_mesh,
      const Direction<dim> direction_to_reconstruct) const;
};

bool operator==(const MonotonisedCentral& /*lhs*/,
                const MonotonisedCentral& /*rhs*/);

bool operator!=(const MonotonisedCentral& lhs, const MonotonisedCentral& rhs);

}  // namespace ForceFree::fd
