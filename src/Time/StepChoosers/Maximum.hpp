// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <utility>

#include "Options/Options.hpp"
#include "Options/String.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/TimeStepRequest.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace StepChoosers {

/// Limits the step size to a constant.
class Maximum : public StepChooser<StepChooserUse::Slab>,
                public StepChooser<StepChooserUse::LtsStep> {
 public:
  /// \cond
  Maximum() = default;
  explicit Maximum(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Maximum);  // NOLINT
  /// \endcond

  static constexpr Options::String help{"Limits the step size to a constant."};

  explicit Maximum(double value);

  using argument_tags = tmpl::list<>;

  std::pair<TimeStepRequest, bool> operator()(double last_step) const;

  bool uses_local_data() const override;
  bool can_be_delayed() const override;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  double value_ = std::numeric_limits<double>::signaling_NaN();
};
}  // namespace StepChoosers

template <>
struct Options::create_from_yaml<StepChoosers::Maximum> {
  template <typename Metavariables>
  static StepChoosers::Maximum create(const Options::Option& options) {
    return create<void>(options);
  }
};

template <>
StepChoosers::Maximum
Options::create_from_yaml<StepChoosers::Maximum>::create<void>(
    const Options::Option& options);
