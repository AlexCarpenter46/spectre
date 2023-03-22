// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>
#include <string>

#include "Options/Options.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class Matrix;
template <size_t Dim>
class Mesh;

/// \endcond

namespace Filters {

/*!
 * \ingroup DiscontinuousGalerkinGroup
 * \brief A cached exponential filter.
 *
 * Applies an exponential filter in each logical direction to each component of
 * the tensors `TagsToFilter`. The exponential filter rescales the 1d modal
 * coefficients \f$c_i\f$ as:
 *
 * \f{align*}{
 *  c_i\to c_i \exp\left[-\alpha_{\mathrm{ef}}
 *   \left(\frac{i}{N}\right)^{2\beta_{\mathrm{ef}}}\right]
 * \f}
 *
 * where \f$N\f$ is the basis degree (number of grid points per element per
 * dimension minus one), \f$\alpha_{\mathrm{ef}}\f$ determines how much the
 * coefficients are rescaled, and \f$\beta_{\mathrm{ef}}\f$ (given by the
 * `HalfPower` option) determines how aggressive/broad the filter is (lower
 * values means filtering more coefficients). Setting
 * \f$\alpha_{\mathrm{ef}}=36\f$ results in effectively zeroing the highest
 * coefficient (in practice it gets rescaled by machine epsilon). The same
 * \f$\alpha_{\mathrm{ef}}\f$ and \f$\beta_{\mathrm{ef}}\f$ are used in each
 * logical direction. For a discussion of filtering see section 5.3 of
 * \cite HesthavenWarburton.
 *
 * #### Design decision:
 *
 * - The reason for the `size_t` template parameter is to allow for different
 * `Alpha` and `HalfPower` parameters for different tensors while still being
 * able to cache the matrices.   If different `Alpha` or `HalfPower` parameters
 * are desired for filtering different tensors, then multiple filters must be
 * inserted into the GlobalCache with different `FilterIndex` values. In
 * the input file these will be specified as `ExpFilterFILTER_INDEX`, e.g.
 * \snippet LinearOperators/Test_Filtering.cpp multiple_exponential_filters
 */
template <size_t FilterIndex>
class Exponential {
 public:
  /// \brief The value of `exp(-alpha)` is what the highest modal coefficient is
  /// rescaled by.
  struct Alpha {
    using type = double;
    static constexpr Options::String help =
        "exp(-alpha) is rescaling of highest coefficient";
    static type lower_bound() { return 0.0; }
  };

  /*!
   * \brief Half of the exponent in the exponential.
   *
   * \f{align*}{
   *  c_i\to c_i \exp\left[-\alpha \left(\frac{i}{N}\right)^{2m}\right]
   * \f}
   */
  struct HalfPower {
    using type = unsigned;
    static constexpr Options::String help =
        "Half of the exponent in the generalized Gaussian";
    static type lower_bound() { return 1; }
  };

  /// \brief Turn the filter off
  struct Enable {
    using type = bool;
    static constexpr Options::String help = {"Enable the filter"};
  };

  using options = tmpl::list<Alpha, HalfPower, Enable>;
  static constexpr Options::String help = {"An exponential filter."};
  static std::string name() {
    return "ExpFilter" + std::to_string(FilterIndex);
  }

  Exponential() = default;

  Exponential(double alpha, unsigned half_power, bool enable);

  /// A cached matrix used to apply the filter to the given mesh
  const Matrix& filter_matrix(const Mesh<1>& mesh) const;

  bool enable() const { return enable_; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 private:
  template <size_t LocalFilterIndex>
  // NOLINTNEXTLINE(readability-redundant-declaration)
  friend bool operator==(const Exponential<LocalFilterIndex>& lhs,
                         const Exponential<LocalFilterIndex>& rhs);

  double alpha_{36.0};
  unsigned half_power_{16};
  bool enable_{true};
};

template <size_t LocalFilterIndex>
bool operator==(const Exponential<LocalFilterIndex>& lhs,
                const Exponential<LocalFilterIndex>& rhs);

template <size_t FilterIndex>
bool operator!=(const Exponential<FilterIndex>& lhs,
                const Exponential<FilterIndex>& rhs);
}  // namespace Filters
