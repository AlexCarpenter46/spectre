// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

// FIXME
#include <cmath>
#include <limits>
#include <pup.h>
#include <random>
#include <utility>

#include "Options/String.hpp"
#include "Time/StepChoosers/StepChooser.hpp"  // IWYU pragma: keep
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
template <size_t VolumeDim>
struct Element;
}  // namespace domain::Tags
/// \endcond

namespace StepChoosers {
/// Changes the step size pseudo-randomly.  Values are distributed
/// uniformly in $\log(dt)$.
template <typename StepChooserUse, size_t VolumeDim>
class Random : public StepChooser<StepChooserUse> {
 public:
  /// \cond
  Random() = default;
  explicit Random(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Random);  // NOLINT
  /// \endcond

  struct Minimum {
    using type = double;
    static constexpr Options::String help{"Minimum value to suggest"};
    static type lower_bound() { return 0.0; }
  };

  struct Maximum {
    using type = double;
    static constexpr Options::String help{"Maximum value to suggest"};
    static type lower_bound() { return 0.0; }
  };

  struct Seed {
    using type = size_t;
    static constexpr Options::String help{"RNG seed"};
  };

  static constexpr Options::String help =
      "Changes the step size pseudo-randomly.";
  using options = tmpl::list<Minimum, Maximum, Seed>;

  explicit Random(const double minimum, const double maximum, const size_t seed,
                  const Options::Context& context = {})
      : minimum_(minimum), maximum_(maximum), seed_(seed) {
    if (minimum_ >= maximum_) {
      PARSE_ERROR(context, "Must have Minimum < Maximum");
    }
  }

  using argument_tags =
      tmpl::list<domain::Tags::Element<VolumeDim>, Tags::TimeStepId>;

  template <typename Element>
  std::pair<double, bool> operator()(
      const Element& element, const TimeStepId& time_step_id,
      const double /*last_step_magnitude*/) const {
    size_t local_seed = seed_;
    boost::hash_combine(local_seed, element.id());
    boost::hash_combine(local_seed, time_step_id);
    std::mt19937_64 rng(local_seed);
    std::uniform_real_distribution<> dist(log(minimum_), log(maximum_));
    for (;;) {
      const double step = exp(dist(rng));
      // Don't produce out-of-range values because of roundoff.
      if (step >= minimum_ and step <= maximum_) {
        return std::make_pair(step, true);
      }
    }
  }

  bool uses_local_data() const override { return false; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    p | minimum_;
    p | maximum_;
    p | seed_;
  }

 private:
  double minimum_ = std::numeric_limits<double>::signaling_NaN();
  double maximum_ = std::numeric_limits<double>::signaling_NaN();
  double seed_ = std::numeric_limits<double>::signaling_NaN();
};

/// \cond
template <typename StepChooserUse, size_t VolumeDim>
PUP::able::PUP_ID Random<StepChooserUse, VolumeDim>::my_PUP_ID = 0;  // NOLINT
/// \endcond
}  // namespace StepChoosers
