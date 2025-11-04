// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Ringdown/StrahlkorperCoefsInRingdownDistortedFrame.hpp"

#include <array>
#include <cstddef>
#include <iostream>
#include <optional>
#include <vector>

#include "DataStructures/Matrix.hpp"
#include "Domain/CoordinateMaps/Distribution.hpp"
#include "Domain/CoordsToDifferentFrame.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Creators/TimeDependentOptions/ExpansionMap.hpp"
#include "Domain/Creators/TimeDependentOptions/FromVolumeFile.hpp"
#include "Domain/Creators/TimeDependentOptions/RotationMap.hpp"
#include "Domain/Creators/TimeDependentOptions/Sphere.hpp"
#include "Domain/Creators/TimeDependentOptions/TranslationMap.hpp"
#include "Domain/StrahlkorperTransformations.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/VolumeData.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/ChangeCenterOfStrahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/IO/ReadSurfaceYlm.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace evolution::Ringdown {
std::pair<std::vector<DataVector>, std::vector<std::array<double, 3>>>
strahlkorper_coefs_in_ringdown_distorted_frame(
    const std::string& path_to_volume_data,
    const std::string& volume_subfile_name,
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name,
    const size_t requested_number_of_times_from_end, const double match_time,
    const double settling_timescale,
    const std::optional<std::array<double, 3>>& exp_func_and_2_derivs,
    const std::optional<std::array<double, 3>>&
        exp_outer_bdry_func_and_2_derivs,
    const std::optional<std::vector<std::array<double, 4>>>&
        rot_func_and_2_derivs,
    const std::optional<std::array<std::array<double, 3>, 3>>&
        trans_func_and_2_derivs) {
  // Read the AhC coefficients from the H5 file
  const std::vector<ylm::Strahlkorper<Frame::Inertial>>& ahc_inertial_h5 =
      ylm::read_surface_ylm<Frame::Inertial>(
          path_to_horizons_h5, surface_subfile_name,
          requested_number_of_times_from_end);
  std::vector<double> ahc_times{};
  {
    // Read the AhC times from the H5 file
    const h5::H5File<h5::AccessType::ReadOnly> ahc_h5_file{path_to_horizons_h5};
    const auto& dat = ahc_h5_file.get<h5::Dat>(surface_subfile_name);
    const Matrix& coefs_for_times = dat.get_data_subset(
        {0}, dat.get_dimensions()[0] - requested_number_of_times_from_end,
        ahc_inertial_h5.size());
    for (size_t i = 0; i < coefs_for_times.rows(); ++i) {
      ahc_times.push_back(coefs_for_times(i, 0));
    }
  }

  const h5::H5File<h5::AccessType::ReadOnly> volume_file{path_to_volume_data};
  const auto& volume_data =
      volume_file.get<h5::VolumeData>(volume_subfile_name);
  const auto obs_ids = volume_data.list_observation_ids();
  size_t obs_id_at_match_time = 0;
  for (const auto obs_id : obs_ids) {
    if (volume_data.get_observation_value(obs_id) == match_time) {
      obs_id_at_match_time = obs_id;
    }
  }

  const auto serialized_inspiral_domain =
      volume_data.get_domain(obs_id_at_match_time);
  if (not serialized_inspiral_domain.has_value()) {
    ERROR("No domain found in volume files at the specified match time.");
  }
  const auto inspiral_domain =
      deserialize<Domain<3>>(serialized_inspiral_domain->data());

  const auto serialized_inspiral_functions_of_time =
      volume_data.get_functions_of_time(obs_id_at_match_time);
  if (not serialized_inspiral_functions_of_time.has_value()) {
    ERROR("No functions of time found in volume files at the match time.");
  }
  const auto inspiral_functions_of_time = deserialize<std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>>(
      serialized_inspiral_functions_of_time->data());

  // Create a time-dependent domain; only the the time-dependent map options
  // matter; the domain is just a spherical shell with inner and outer
  // radii chosen so any conceivable common horizon will fit between them.
  const auto expansion_map_options =
      exp_func_and_2_derivs.has_value()
          ? domain::creators::time_dependent_options::ExpansionMapOptions<
                true>{exp_func_and_2_derivs.value(), settling_timescale,
                      exp_outer_bdry_func_and_2_derivs.value(),
                      settling_timescale}
          : std::optional<domain::creators::time_dependent_options::
                              ExpansionMapOptions<true>>{};
  const auto rotation_map_options =
      rot_func_and_2_derivs.has_value()
          ? domain::creators::time_dependent_options::RotationMapOptions<
                true>{rot_func_and_2_derivs.value(), settling_timescale}
          : std::optional<domain::creators::time_dependent_options::
                              RotationMapOptions<true>>{};

  const auto& translation_fot_from_volume =
      trans_func_and_2_derivs.has_value()
          ? domain::creators::time_dependent_options::FromVolumeFile(
                path_to_volume_data, volume_subfile_name, true)
          : std::optional<
                domain::creators::time_dependent_options::FromVolumeFile>{};
  const domain::creators::sphere::TimeDependentMapOptions
      time_dependent_map_options{match_time,
                                 std::nullopt,
                                 rotation_map_options,
                                 expansion_map_options,
                                 translation_fot_from_volume,
                                 true};
  const domain::creators::Sphere domain_creator{
      0.01,
      200.0,
      // nullptr because no boundary condition
      domain::creators::Sphere::Excision{nullptr},
      static_cast<size_t>(0),
      static_cast<size_t>(5),
      false,
      std::nullopt,
      {50.0},
      domain::CoordinateMaps::Distribution::Linear,
      ShellWedges::All,
      time_dependent_map_options};
  const auto temporary_domain = domain_creator.create_domain();
  const auto functions_of_time = domain_creator.functions_of_time();
  // Loop over the selected horizons, transforming each to the
  // ringdown distorted frame
  std::vector<DataVector> ahc_ringdown_distorted_coefs{};
  std::vector<std::array<double, 3>> ahc_inertial_centers{};
  // Here we transform the inertial strahlkorper into the ringdown distorted
  // frame. In order to do this, the inertial coords of the strahlkorper are
  // mapped to the logical frame to determine which block map to use, and then
  // into the ringdown distorted frame. This technically requires a shape map,
  // however, at this time we do not yet know the shape map for the ringdown
  // domain so we can use an identity instead because we are only concerned with
  // the correct ringdown distorted frame, not the correct ringdown grid frame.
  // To avoid an unnecessary identity shape map, we omit it in the domain above.
  // This now makes the grid frame of this temporary domain equivalent to the
  // true ringdown distorted frame we are after. This is why we map the
  // strahlkorper into the "grid" frame instead of the "distorted" frame. It is
  // a simplification to avoid an unnecessary identity shape map.
  ylm::Strahlkorper<Frame::Grid> distorted_ahc;
  ylm::Strahlkorper<Frame::Inertial> inertial_ahc;
  for (size_t i = 0; i < requested_number_of_times_from_end; ++i) {
    if (gsl::at(ahc_times, i) <= match_time) {
      std::cout << "Time: " << ahc_times[i] << '\n';
      std::cout << "Original inertial frame strahlkorper physical center from file: " <<
      gsl::at(ahc_inertial_h5, i).physical_center() << '\n';
      std::cout << "Original inertial frame strahlkorper expansion center from file: " <<
      gsl::at(ahc_inertial_h5, i).expansion_center() << '\n';
      strahlkorper_in_different_frame(
          make_not_null(&distorted_ahc), gsl::at(ahc_inertial_h5, i),
          temporary_domain, functions_of_time, gsl::at(ahc_times, i));
      // Relative tolerance is set to the value used in SpEC
      ylm::change_expansion_center_of_strahlkorper_to_physical(
          make_not_null(&distorted_ahc), 1e-7);
      ahc_ringdown_distorted_coefs.push_back(distorted_ahc.coefficients());
      strahlkorper_in_different_frame(
          make_not_null(&inertial_ahc), distorted_ahc, inspiral_domain,
          inspiral_functions_of_time, gsl::at(ahc_times, i));
      std::cout
          << "Inertial center point from strahl in diff frame phys center: "
          << inertial_ahc.physical_center() << '\n';
      std::cout
          << "Inertial center point from strahl in diff frame exp center: "
          << inertial_ahc.expansion_center() << '\n';
      tnsr::I<DataVector, 3, ::Frame::Grid> grid_center_point{
          DataVector{1, 0.0}};
      grid_center_point[0] = distorted_ahc.expansion_center()[0];
      grid_center_point[1] = distorted_ahc.expansion_center()[1];
      grid_center_point[2] = distorted_ahc.expansion_center()[2];
      std::cout << "Inspiral inertial center point: "
                << gsl::at(ahc_inertial_h5, i).physical_center() << '\n';
      std::cout << "grid center point: " << grid_center_point << '\n';
      tnsr::I<DataVector, 3, ::Frame::Inertial> inertial_center_point{
          DataVector{1, 0.0}};
      // The center point is mapped back to the inspiral inertial frame so that
      // the center of AhC is the same at the match time.
      coords_to_different_frame(
          make_not_null(&inertial_center_point), grid_center_point,
          inspiral_domain, inspiral_functions_of_time, gsl::at(ahc_times, i));
      std::cout << "Inspiral inertial center point from coords to diff frame: "
                << inertial_center_point << '\n';
      // This minus sign was found by trial and error.
    //   ahc_inertial_centers.push_back(
    //       std::array<double, 3>{-1.0 * get<0>(inertial_center_point)[0],
    //                             -1.0 * get<1>(inertial_center_point)[0],
    //                             -1.0 * get<2>(inertial_center_point)[0]});
      ahc_inertial_centers.push_back(gsl::at(ahc_inertial_h5, i).physical_center());
    }
  }

  return std::pair{ahc_ringdown_distorted_coefs, ahc_inertial_centers};
}
}  // namespace evolution::Ringdown
