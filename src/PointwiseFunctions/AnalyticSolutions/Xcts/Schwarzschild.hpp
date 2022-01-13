// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <ostream>

#include "DataStructures/CachedTempBuffer.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/CommonVariables.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags/Conformal.hpp"
#include "PointwiseFunctions/InitialDataUtilities/AnalyticSolution.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace Xcts::Solutions {

/// Various coordinate systems in which to express the Schwarzschild solution
enum class SchwarzschildCoordinates {
  /*!
   * \brief Isotropic Schwarzschild coordinates
   *
   * These arise from the canonical Schwarzschild coordinates by the radial
   * transformation
   *
   * \f{equation}
   * r = \bar{r}\left(1+\frac{M}{2\bar{r}}\right)^2
   * \f}
   *
   * (Eq. (1.61) in \cite BaumgarteShapiro) where \f$r\f$ is the canonical
   * Schwarzschild radius, also referred to as "areal" radius because it is
   * defined such that spheres with constant \f$r\f$ have the area \f$4\pi
   * r^2\f$, and \f$\bar{r}\f$ is the "isotropic" radius. In the isotropic
   * radius the Schwarzschild spatial metric is conformally flat:
   *
   * \f{equation}
   * \gamma_{ij}=\psi^4\eta_{ij} \quad \text{with conformal factor} \quad
   * \psi=1+\frac{M}{2\bar{r}}
   * \f}
   *
   * (Table 2.1 in \cite BaumgarteShapiro). Its lapse transforms to
   *
   * \f{equation}
   * \alpha=\frac{1-M/(2\bar{r})}{1+M/(2\bar{r})}
   * \f}
   *
   * and the shift vanishes (\f$\beta^i=0\f$) as it does in areal Schwarzschild
   * coordinates. The solution also remains maximally sliced, i.e. \f$K=0\f$.
   *
   * The Schwarzschild horizon in these coordinates is at
   * \f$\bar{r}=\frac{M}{2}\f$ due to the radial transformation from \f$r=2M\f$.
   */
  Isotropic,
};

std::ostream& operator<<(std::ostream& os, SchwarzschildCoordinates coords);

}  // namespace Xcts::Solutions

template <>
struct Options::create_from_yaml<Xcts::Solutions::SchwarzschildCoordinates> {
  template <typename Metavariables>
  static Xcts::Solutions::SchwarzschildCoordinates create(
      const Options::Option& options) {
    return create<void>(options);
  }
};

template <>
Xcts::Solutions::SchwarzschildCoordinates
Options::create_from_yaml<Xcts::Solutions::SchwarzschildCoordinates>::create<
    void>(const Options::Option& options);

namespace Xcts::Solutions {

namespace detail {

struct SchwarzschildImpl {
  struct Mass {
    using type = double;
    static constexpr Options::String help = "Mass parameter M";
  };

  struct CoordinateSystem {
    static std::string name() { return "Coordinates"; }
    using type = SchwarzschildCoordinates;
    static constexpr Options::String help =
        "The coordinate system used to describe the solution";
  };

  using options = tmpl::list<Mass, CoordinateSystem>;
  static constexpr Options::String help{
      "Schwarzschild spacetime in general relativity"};

  SchwarzschildImpl() = default;
  SchwarzschildImpl(const SchwarzschildImpl&) = default;
  SchwarzschildImpl& operator=(const SchwarzschildImpl&) = default;
  SchwarzschildImpl(SchwarzschildImpl&&) = default;
  SchwarzschildImpl& operator=(SchwarzschildImpl&&) = default;
  ~SchwarzschildImpl() = default;

  explicit SchwarzschildImpl(double mass,
                             SchwarzschildCoordinates coordinate_system);

  /// The mass parameter \f$M\f$.
  double mass() const;

  SchwarzschildCoordinates coordinate_system() const;

  /// The radius of the Schwarzschild horizon in the given coordinates.
  double radius_at_horizon() const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 protected:
  double mass_{std::numeric_limits<double>::signaling_NaN()};
  SchwarzschildCoordinates coordinate_system_{};
};

bool operator==(const SchwarzschildImpl& lhs, const SchwarzschildImpl& rhs);

bool operator!=(const SchwarzschildImpl& lhs, const SchwarzschildImpl& rhs);

template <typename DataType>
using SchwarzschildVariablesCache = cached_temp_buffer_from_typelist<
    tmpl::push_front<
        common_tags<DataType>,
        gr::Tags::Conformal<gr::Tags::EnergyDensity<DataType>, 0>,
        gr::Tags::Conformal<gr::Tags::StressTrace<DataType>, 0>,
        gr::Tags::Conformal<
            gr::Tags::MomentumDensity<3, Frame::Inertial, DataType>, 0>>>;

template <typename DataType>
struct SchwarzschildVariables
    : CommonVariables<DataType, SchwarzschildVariablesCache<DataType>> {
  static constexpr size_t Dim = 3;
  static constexpr int ConformalMatterScale = 0;
  using Cache = SchwarzschildVariablesCache<DataType>;
  using Base = CommonVariables<DataType, SchwarzschildVariablesCache<DataType>>;
  using Base::operator();

  SchwarzschildVariables(
      std::optional<std::reference_wrapper<const Mesh<Dim>>> local_mesh,
      std::optional<std::reference_wrapper<const InverseJacobian<
          DataType, Dim, Frame::ElementLogical, Frame::Inertial>>>
          local_inv_jacobian,
      const tnsr::I<DataType, 3>& local_x, const double local_mass,
      const SchwarzschildCoordinates local_coordinate_system)
      : Base(std::move(local_mesh), std::move(local_inv_jacobian)),
        x(local_x),
        mass(local_mass),
        coordinate_system(local_coordinate_system) {}

  const tnsr::I<DataType, 3>& x;
  double mass;
  SchwarzschildCoordinates coordinate_system;

  void operator()(gsl::not_null<tnsr::ii<DataType, 3>*> conformal_metric,
                  gsl::not_null<Cache*> cache,
                  Tags::ConformalMetric<DataType, 3, Frame::Inertial> /*meta*/)
      const override;
  void operator()(
      gsl::not_null<tnsr::II<DataType, 3>*> inv_conformal_metric,
      gsl::not_null<Cache*> cache,
      Tags::InverseConformalMetric<DataType, 3, Frame::Inertial> /*meta*/)
      const override;
  void operator()(
      gsl::not_null<tnsr::ijj<DataType, 3>*> deriv_conformal_metric,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::ConformalMetric<DataType, 3, Frame::Inertial>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;
  void operator()(
      gsl::not_null<Scalar<DataType>*> trace_extrinsic_curvature,
      gsl::not_null<Cache*> cache,
      gr::Tags::TraceExtrinsicCurvature<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::i<DataType, 3>*> trace_extrinsic_curvature_gradient,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataType>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;
  void operator()(
      gsl::not_null<Scalar<DataType>*> dt_trace_extrinsic_curvature,
      gsl::not_null<Cache*> cache,
      ::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataType>> /*meta*/)
      const override;
  void operator()(gsl::not_null<Scalar<DataType>*> conformal_factor,
                  gsl::not_null<Cache*> cache,
                  Tags::ConformalFactor<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::i<DataType, 3>*> conformal_factor_gradient,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Xcts::Tags::ConformalFactor<DataType>, tmpl::size_t<3>,
                    Frame::Inertial> /*meta*/) const override;
  void operator()(gsl::not_null<Scalar<DataType>*> lapse,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Lapse<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor,
      gsl::not_null<Cache*> cache,
      Tags::LapseTimesConformalFactor<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::i<DataType, 3>*>
          lapse_times_conformal_factor_gradient,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::LapseTimesConformalFactor<DataType>, tmpl::size_t<3>,
                    Frame::Inertial> /*meta*/) const override;
  void operator()(gsl::not_null<tnsr::I<DataType, 3>*> shift_background,
                  gsl::not_null<Cache*> cache,
                  Tags::ShiftBackground<DataType, 3, Frame::Inertial> /*meta*/)
      const override;
  void operator()(gsl::not_null<tnsr::II<DataType, 3, Frame::Inertial>*>
                      longitudinal_shift_background_minus_dt_conformal_metric,
                  gsl::not_null<Cache*> cache,
                  Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
                      DataType, 3, Frame::Inertial> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::I<DataType, 3>*> shift_excess,
      gsl::not_null<Cache*> cache,
      Tags::ShiftExcess<DataType, 3, Frame::Inertial> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::ii<DataType, 3>*> shift_strain,
      gsl::not_null<Cache*> cache,
      Tags::ShiftStrain<DataType, 3, Frame::Inertial> /*meta*/) const override;
  void operator()(gsl::not_null<Scalar<DataType>*> energy_density,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Conformal<gr::Tags::EnergyDensity<DataType>,
                                      ConformalMatterScale> /*meta*/) const;
  void operator()(gsl::not_null<Scalar<DataType>*> stress_trace,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Conformal<gr::Tags::StressTrace<DataType>,
                                      ConformalMatterScale> /*meta*/) const;
  void operator()(gsl::not_null<tnsr::I<DataType, 3>*> momentum_density,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Conformal<
                      gr::Tags::MomentumDensity<3, Frame::Inertial, DataType>,
                      ConformalMatterScale> /*meta*/) const;
};

}  // namespace detail

/*!
 * \brief Schwarzschild spacetime in general relativity
 *
 * This class implements the Schwarzschild solution with mass parameter
 * \f$M\f$ in various coordinate systems. See the entries of the
 * `Xcts::Solutions::SchwarzschildCoordinates` enum for the available coordinate
 * systems and for the solution variables in the respective coordinates.
 */
class Schwarzschild : public elliptic::analytic_data::AnalyticSolution,
                      public detail::SchwarzschildImpl {
 public:
  Schwarzschild() = default;
  Schwarzschild(const Schwarzschild&) = default;
  Schwarzschild& operator=(const Schwarzschild&) = default;
  Schwarzschild(Schwarzschild&&) = default;
  Schwarzschild& operator=(Schwarzschild&&) = default;
  ~Schwarzschild() = default;

  using SchwarzschildImpl::SchwarzschildImpl;

  /// \cond
  explicit Schwarzschild(CkMigrateMessage* m)
      : elliptic::analytic_data::AnalyticSolution(m) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Schwarzschild);
  /// \endcond

  template <typename DataType, typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables(
      const tnsr::I<DataType, 3, Frame::Inertial>& x,
      tmpl::list<RequestedTags...> /*meta*/) const {
    using VarsComputer = detail::SchwarzschildVariables<DataType>;
    typename VarsComputer::Cache cache{get_size(*x.begin())};
    const VarsComputer computer{std::nullopt, std::nullopt, x, mass_,
                                coordinate_system_};
    return {cache.get_var(computer, RequestedTags{})...};
  }

  template <typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables(
      const tnsr::I<DataVector, 3, Frame::Inertial>& x, const Mesh<3>& mesh,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& inv_jacobian,
      tmpl::list<RequestedTags...> /*meta*/) const {
    using VarsComputer = detail::SchwarzschildVariables<DataVector>;
    typename VarsComputer::Cache cache{get_size(*x.begin())};
    const VarsComputer computer{mesh, inv_jacobian, x, mass_,
                                coordinate_system_};
    return {cache.get_var(computer, RequestedTags{})...};
  }

  void pup(PUP::er& p) override {
    elliptic::analytic_data::AnalyticSolution::pup(p);
    detail::SchwarzschildImpl::pup(p);
  }
};

}  // namespace Xcts::Solutions
