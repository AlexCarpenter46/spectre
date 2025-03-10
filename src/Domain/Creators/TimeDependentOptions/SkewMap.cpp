// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependentOptions/SkewMap.hpp"

#include <array>
#include <string>
#include <variant>

#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/TimeDependentOptions/FromVolumeFile.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"

namespace domain::creators::time_dependent_options {
template <bool AllowReplay>
SkewMapOptions<AllowReplay>::SkewMapOptions(
    const std::array<double, 3>& initial_angles_y_in,
    const std::array<double, 3>& initial_angles_z_in)
    : initial_angles_y(initial_angles_y_in),
      initial_angles_z(initial_angles_z_in) {}

template <bool AllowReplay>
std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> get_skew(
    const std::variant<SkewMapOptions<AllowReplay>,
                       FromVolumeFile<AllowReplay>>& skew_map_options,
    const double initial_time, const double expiration_time) {
  const std::string name{"Skew"};
  std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> result{};

  if (std::holds_alternative<FromVolumeFile<AllowReplay>>(skew_map_options)) {
    const auto& from_vol_file =
        std::get<FromVolumeFile<AllowReplay>>(skew_map_options);
    auto volume_fot =
        from_vol_file.retrieve_function_of_time({name}, initial_time);

    // It must be a PiecewisePolynomial
    if (UNLIKELY(dynamic_cast<domain::FunctionsOfTime::PiecewisePolynomial<2>*>(
                     volume_fot.at(name).get()) == nullptr)) {
      ERROR_NO_TRACE(
          "Skew function of time read from volume data is not a "
          "PiecewisePolynomial<2>. Cannot use it to initialize the skew map.");
    }

    if (from_vol_file.replay()) {
      result = std::move(volume_fot.at(name));
    } else {
      result =
          volume_fot.at(name)->create_at_time(initial_time, expiration_time);
    }
  } else if (std::holds_alternative<SkewMapOptions<AllowReplay>>(
                 skew_map_options)) {
    const auto& hard_coded_options =
        std::get<SkewMapOptions<AllowReplay>>(skew_map_options);

    result = std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
        initial_time,
        std::array{
            DataVector{hard_coded_options.initial_angles_y[0],
                       hard_coded_options.initial_angles_z[0]},
            DataVector{hard_coded_options.initial_angles_y[1],
                       hard_coded_options.initial_angles_z[1]},
            DataVector{hard_coded_options.initial_angles_y[2],
                       hard_coded_options.initial_angles_z[2]},
        },
        expiration_time);
  } else {
    ERROR("Unknown SkewMap.");
  }

  return result;
}

template struct SkewMapOptions<true>;
template struct SkewMapOptions<false>;
template std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> get_skew(
    const std::variant<SkewMapOptions<true>, FromVolumeFile<true>>&
        skew_map_options,
    const double initial_time, const double expiration_time);
template std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> get_skew(
    const std::variant<SkewMapOptions<false>, FromVolumeFile<false>>&
        skew_map_options,
    const double initial_time, const double expiration_time);
}  // namespace domain::creators::time_dependent_options
