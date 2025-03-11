// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Ringdown/MinimumAhCExcisionRadius.hpp"

#include <array>
#include <atomic>
#include <cstddef>
#include <optional>
#include <vector>

#include <iostream>

#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/CoordinateMaps/Distribution.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Domain/CoordsToDifferentFrame.hpp"
#include "Domain/Creators/BinaryCompactObject.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Creators/TimeDependentOptions/BinaryCompactObject.hpp"
#include "Domain/Creators/TimeDependentOptions/ExpansionMap.hpp"
#include "Domain/Creators/TimeDependentOptions/RotationMap.hpp"
#include "Domain/Creators/TimeDependentOptions/ShapeMap.hpp"
#include "Domain/Creators/TimeDependentOptions/Sphere.hpp"
#include "Domain/StrahlkorperTransformations.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/VolumeData.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/IO/ReadSurfaceYlm.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/StrahlkorperFunctions.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace evolution::Ringdown {
double minimum_ahc_excision_radius(
    const std::string& path_to_volume_data,
    const std::string& volume_subfile_name,
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name,
    const std::string& path_to_AhC_distorted_h5,
    const std::vector<std::string>& AhC_distorted_subfile_names,
    size_t requested_number_of_times_from_end, double match_time,
    double settling_timescale, double excision_A_radius,
    double excision_B_radius, std::array<double, 3> excision_A_center,
    std::array<double, 3> excision_B_center,
    const std::optional<std::array<double, 3>>& exp_func_and_2_derivs,
    const std::optional<std::array<double, 3>>&
        exp_outer_bdry_func_and_2_derivs,
    const std::optional<std::vector<std::array<double, 4>>>&
        rot_func_and_2_derivs,
    const std::optional<std::array<std::array<double, 3>, 3>>&
        trans_func_and_2_derivs) {
  // Read the AhC coefficients from the H5 file
  const std::vector<ylm::Strahlkorper<Frame::Inertial>>& ahc_inertial =
      ylm::read_surface_ylm<Frame::Inertial>(
          path_to_horizons_h5, surface_subfile_name,
          requested_number_of_times_from_end);

  ylm::Strahlkorper<Frame::Inertial> ahc_inertial_at_match_time{};
  std::vector<double> ahc_times{};
  // Read the AhC times from the H5 file
  const h5::H5File<h5::AccessType::ReadOnly> ahc_h5_file{path_to_horizons_h5};
  const auto& dat = ahc_h5_file.get<h5::Dat>(surface_subfile_name);
  const Matrix& coefs_for_times = dat.get_data_subset(
      {0}, dat.get_dimensions()[0] - requested_number_of_times_from_end,
      ahc_inertial.size());
  for (size_t i = 0; i < coefs_for_times.rows(); ++i) {
    ahc_times.push_back(coefs_for_times(i, 0));
  }
  for (size_t i = 0; i < ahc_times.size(); i++) {
    if (gsl::at(ahc_times, i) == match_time) {
      ahc_inertial_at_match_time = gsl::at(ahc_inertial, i);
    }
  }

  // Create a time-dependent domain; only the the time-dependent map options
  // matter; the domain is just a spherical shell with inner and outer
  // radii chosen so any conceivable common horizon will fit between them.

  // Make sure to change l_max at some point :) WhY 20??? Maybe make parameter
  // but check pl script first.
  const auto shape_map_options =
      domain::creators::time_dependent_options::ShapeMapOptions<
          false, domain::ObjectLabel::None>{
          33,
          domain::creators::time_dependent_options::YlmsFromFile{
              path_to_AhC_distorted_h5, AhC_distorted_subfile_names, match_time,
              1.e-10, true, true},
          std::array<double, 3>{0.0, -1.0, 0.0}};

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
  const auto translation_map_options =
      trans_func_and_2_derivs.has_value()
          ? domain::creators::sphere::TimeDependentMapOptions::
                TranslationMapOptions{trans_func_and_2_derivs.value()}
          : std::optional<domain::creators::sphere::TimeDependentMapOptions::
                              TranslationMapOptions>{};

  const domain::creators::sphere::TimeDependentMapOptions
      ringdown_time_dependent_map_options{match_time,
                                          shape_map_options,
                                          rotation_map_options,
                                          expansion_map_options,
                                          translation_map_options,
                                          true};

  const double rAH = ahc_inertial_at_match_time.average_radius();

  std::cout << "rAH: " << rAH << std::endl;
  std::cout << "match time: " << match_time << std::endl;
  std::cout << "excision A radius: " << excision_A_radius << std::endl;
  std::cout << "excision B radius: " << excision_B_radius << std::endl;
  std::cout << "excision A center: " << excision_A_center << std::endl;
  std::cout << "excision B center: " << excision_B_center << std::endl;

  std::vector<double> vol_times{};
  // Read the AhC times from the H5 file
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
    ERROR("No domain in volume files. Goodnight");
  }
  const auto inspiral_domain =
      deserialize<Domain<3>>(serialized_inspiral_domain->data());

  const auto serialized_inspiral_functions_of_time =
      volume_data.get_functions_of_time(obs_id_at_match_time);
  if (not serialized_inspiral_functions_of_time.has_value()) {
    ERROR("No functions of time in volume files. Goodnight");
  }
  const auto inspiral_functions_of_time = deserialize<std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>>(
      serialized_inspiral_functions_of_time->data());

  // Poor-man's mass ratio.
  double q = abs(excision_B_center[0] / excision_A_center[0]);
  std::cout << "Poor man's mass ratio: " << q << std::endl;
  double eps = 1e-3 / q;
  bool outer_converged = false;
  size_t current_outer_iteration = 0;
  size_t max_iterations = 10;
  size_t current_l_max = 20;
  double rminfac = 0.94;

  // Start of outer loop
  while (not outer_converged and current_outer_iteration < max_iterations) {
    // Section for constructing strahlkorpers for AhA/AhB
    ylm::Strahlkorper<Frame::Grid> excision_a_inspiral_grid(
        current_l_max, excision_A_radius, excision_A_center);
    ylm::Strahlkorper<Frame::Grid> excision_b_inspiral_grid(
        current_l_max, excision_B_radius, excision_B_center);

    tnsr::I<DataVector, 3, Frame::Inertial> excision_a_inspiral_inertial_points{
        get<0>(ylm::cartesian_coords(excision_a_inspiral_grid)).size()};
    tnsr::I<DataVector, 3, Frame::Inertial> excision_b_inspiral_inertial_points{
        get<0>(ylm::cartesian_coords(excision_b_inspiral_grid)).size()};
    coords_to_different_frame(
        make_not_null(&excision_a_inspiral_inertial_points),
        ylm::cartesian_coords(excision_a_inspiral_grid), inspiral_domain,
        inspiral_functions_of_time, match_time);
    coords_to_different_frame(
        make_not_null(&excision_b_inspiral_inertial_points),
        ylm::cartesian_coords(excision_b_inspiral_grid), inspiral_domain,
        inspiral_functions_of_time, match_time);

    // Loop section finding the correct rmin_fac
    double old_rminfac = rminfac;
    size_t current_inner_iteration = 0;
    bool inner_converged = false;

    // Start of the inner loop!
    while (not inner_converged and current_inner_iteration < max_iterations) {
      const domain::creators::Sphere ringdown_rmin_domain_creator{
          rAH * rminfac,
          200.0,
          // nullptr because no boundary condition
          domain::creators::Sphere::Excision{nullptr},
          static_cast<size_t>(0),
          static_cast<size_t>(5),
          false,
          std::nullopt,
          // make this an option at some point soon :)
          {50.0},
          domain::CoordinateMaps::Distribution::Linear,
          ShellWedges::All,
          ringdown_time_dependent_map_options};

      const auto temporary_ringdown_rmin_domain =
          ringdown_rmin_domain_creator.create_domain();
      const auto ringdown_rmin_functions_of_time =
          ringdown_rmin_domain_creator.functions_of_time();
      const auto exc_a_block_logical = block_logical_coordinates(
          temporary_ringdown_rmin_domain, excision_a_inspiral_inertial_points,
          match_time, ringdown_rmin_functions_of_time);
      const auto exc_b_block_logical = block_logical_coordinates(
          temporary_ringdown_rmin_domain, excision_b_inspiral_inertial_points,
          match_time, ringdown_rmin_functions_of_time);

      tnsr::I<DataVector, 3, Frame::Inertial> excision_a_ringdown_grid_points{
          get<0>(ylm::cartesian_coords(excision_a_inspiral_grid)).size()};
      tnsr::I<DataVector, 3, Frame::Inertial> excision_b_ringdown_grid_points{
          get<0>(ylm::cartesian_coords(excision_b_inspiral_grid)).size()};

      double min_excision_radius = 0.0;
      tnsr::I<double, 3, Frame::Inertial> exc_a_inspiral_inertial_point{};
      tnsr::I<double, 3, Frame::Grid> exc_a_ringdown_grid_point{};
      tnsr::I<double, 3, Frame::Inertial> exc_b_inspiral_inertial_point{};
      tnsr::I<double, 3, Frame::Grid> exc_b_ringdown_grid_point{};
      for (size_t s = 0; s < get<0>(excision_a_inspiral_inertial_points).size();
           ++s) {
        get<0>(exc_a_inspiral_inertial_point) =
            get<0>(excision_a_inspiral_inertial_points)[s];
        get<1>(exc_a_inspiral_inertial_point) =
            get<1>(excision_a_inspiral_inertial_points)[s];
        get<2>(exc_a_inspiral_inertial_point) =
            get<2>(excision_a_inspiral_inertial_points)[s];
        get<0>(exc_b_inspiral_inertial_point) =
            get<0>(excision_b_inspiral_inertial_points)[s];
        get<1>(exc_b_inspiral_inertial_point) =
            get<1>(excision_b_inspiral_inertial_points)[s];
        get<2>(exc_b_inspiral_inertial_point) =
            get<2>(excision_b_inspiral_inertial_points)[s];

        if (exc_a_block_logical[s].has_value()) {
          const auto& block_id_and_coords = exc_a_block_logical[s].value();
          const auto& block = temporary_ringdown_rmin_domain
                                  .blocks()[block_id_and_coords.id.get_index()];
          const auto& grid_to_inertial_map =
              block.moving_mesh_grid_to_inertial_map();
          const auto inv_point = grid_to_inertial_map.inverse(
              exc_a_inspiral_inertial_point, match_time,
              ringdown_rmin_functions_of_time);
          if (inv_point.has_value()) {
            exc_a_ringdown_grid_point = inv_point.value();
          } else {
            ERROR("Map from Frame::Inertial to Frame::Grid is not invertible");
          }
        } else {
          // This point is inside the excision which is why it couldn't be
          // mapped to a block. Now we must use the shape map inverse to
          // determine where this point is.
          const auto inv_point =
              temporary_ringdown_rmin_domain.excision_spheres()
                  .at("ExcisionSphere")
                  .moving_mesh_grid_to_inertial_map()
                  .inverse(exc_a_inspiral_inertial_point, match_time,
                           ringdown_rmin_functions_of_time);
          if (inv_point.has_value()) {
            get<0>(exc_a_ringdown_grid_point) = inv_point.value()[0];
            get<1>(exc_a_ringdown_grid_point) = inv_point.value()[1];
            get<2>(exc_a_ringdown_grid_point) = inv_point.value()[2];
          } else {
            ERROR("Map from Frame::Inertial to Frame::Grid is not invertible");
          }
        }
        if (exc_b_block_logical[s].has_value()) {
          const auto& block_id_and_coords = exc_b_block_logical[s].value();
          const auto& block = temporary_ringdown_rmin_domain
                                  .blocks()[block_id_and_coords.id.get_index()];
          const auto& grid_to_inertial_map =
              block.moving_mesh_grid_to_inertial_map();
          const auto inv_point = grid_to_inertial_map.inverse(
              exc_b_inspiral_inertial_point, match_time,
              ringdown_rmin_functions_of_time);
          if (inv_point.has_value()) {
            exc_b_ringdown_grid_point = inv_point.value();
          } else {
            ERROR("Map from Frame::Inertial to Frame::Grid is not invertible");
          }
        } else {
          // This point is inside the excision which is why it couldn't be
          // mapped to a block. Now we must use the shape map inverse to
          // determine where this point is.
          const auto inv_point =
              temporary_ringdown_rmin_domain.excision_spheres()
                  .at("ExcisionSphere")
                  .moving_mesh_grid_to_inertial_map()
                  .inverse(exc_b_inspiral_inertial_point, match_time,
                           ringdown_rmin_functions_of_time);
          if (inv_point.has_value()) {
            get<0>(exc_b_ringdown_grid_point) = inv_point.value()[0];
            get<1>(exc_b_ringdown_grid_point) = inv_point.value()[1];
            get<2>(exc_b_ringdown_grid_point) = inv_point.value()[2];
          } else {
            ERROR("Map from Frame::Inertial to Frame::Grid is not invertible");
          }
        }
        const double excision_a_point_radius =
            sqrt(square(get<0>(exc_a_ringdown_grid_point) -
                        trans_func_and_2_derivs.value()[0][0]) +
                 square(get<1>(exc_a_ringdown_grid_point) -
                        trans_func_and_2_derivs.value()[0][1]) +
                 square(get<2>(exc_a_ringdown_grid_point) -
                        trans_func_and_2_derivs.value()[0][2]));
        const double excision_b_point_radius =
            sqrt(square(get<0>(exc_b_ringdown_grid_point) -
                        trans_func_and_2_derivs.value()[0][0]) +
                 square(get<1>(exc_b_ringdown_grid_point) -
                        trans_func_and_2_derivs.value()[0][1]) +
                 square(get<2>(exc_b_ringdown_grid_point) -
                        trans_func_and_2_derivs.value()[0][2]));

        if (excision_a_point_radius > min_excision_radius) {
          min_excision_radius = excision_a_point_radius;
        }
        if (excision_b_point_radius > min_excision_radius) {
          min_excision_radius = excision_b_point_radius;
        }
      }
      const double min_rminfac = min_excision_radius / rAH;
      rminfac = 1.0 - 0.25 * (1.0 - min_rminfac);
      if (current_inner_iteration != 0 and
          abs(rminfac - old_rminfac) <= 0.5 * eps) {
        inner_converged = true;
      }
        current_inner_iteration++;
        old_rminfac = rminfac;
    }
    if (current_outer_iteration != 0 and
        abs(rminfac - old_rminfac) <= 0.5 * eps) {
      outer_converged = true;
    }
      current_outer_iteration++;
      // Increment l max by 6 every iteration.
      current_l_max += 6;
      if (current_outer_iteration > max_iterations) {
        ERROR(
            "Max Iterations for finding a suitable excision radius exceeded. "
            "Going to sleep.");
      }
  }
  double safe_rminfac = rminfac / eps;
  safe_rminfac *= eps;
  safe_rminfac += eps;
  if (safe_rminfac - rminfac < 0.5 * eps) {
    safe_rminfac += eps;
  }
  const double excision_radius =
      rminfac > safe_rminfac ? rAH * rminfac : rAH * safe_rminfac;
  std::cout << "AhC Excision Radius: " << excision_radius << std::endl;

  return excision_radius;
}
}  // namespace evolution::Ringdown
