// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/CoordinateMaps/CoordinateMap.hpp"

#include <algorithm>
#include <array>
#include <boost/optional.hpp>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <pup.h>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/Tensor/Identity.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Parallel/CharmPupable.hpp"
#include "Parallel/PupStlCpp11.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/Tuple.hpp"

/// \cond
namespace domain {
// define type-trait to check for time-dependent mapping
namespace CoordinateMap_detail {
template <typename T>
using is_map_time_dependent_t = tt::is_callable_t<
    T, std::array<std::decay_t<T>, std::decay_t<T>::dim>, double,
    std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>>;

template <typename T, size_t Dim, typename Map>
void apply_map(
    const gsl::not_null<std::array<T, Dim>*> t_map_point, const Map& the_map,
    const double /*t*/,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
    /*functions_of_time*/,
    const std::false_type /*is_time_independent*/) {
  if (LIKELY(not the_map.is_identity())) {
    *t_map_point = the_map(*t_map_point);
  }
}

template <typename T, size_t Dim, typename Map>
void apply_map(
    const gsl::not_null<std::array<T, Dim>*> t_map_point, const Map& the_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time,
    const std::true_type
    /*is_time_dependent*/) {
  *t_map_point = the_map(*t_map_point, t, functions_of_time);
}
}  // namespace CoordinateMap_detail

template <typename SourceFrame, typename TargetFrame, typename... Maps>
CoordinateMap<SourceFrame, TargetFrame, Maps...>::CoordinateMap(Maps... maps)
    : maps_(std::move(maps)...) {}

template <typename SourceFrame, typename TargetFrame, typename... Maps>
template <typename T, size_t... Is>
tnsr::I<T, CoordinateMap<SourceFrame, TargetFrame, Maps...>::dim, TargetFrame>
CoordinateMap<SourceFrame, TargetFrame, Maps...>::call_impl(
    tnsr::I<T, dim, SourceFrame>&& source_point, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time,
    std::index_sequence<Is...> /*meta*/) const noexcept {
  std::array<T, dim> mapped_point = make_array<T, dim>(std::move(source_point));

  EXPAND_PACK_LEFT_TO_RIGHT(make_overloader(
      [](const auto& the_map, std::array<T, dim>& point, const double /*t*/,
         const std::unordered_map<
             std::string,
             std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
         /*funcs_of_time*/,
         const std::false_type /*is_time_independent*/) noexcept {
        if (LIKELY(not the_map.is_identity())) {
          point = the_map(point);
        }
      },
      [](const auto& the_map, std::array<T, dim>& point, const double t,
         const std::unordered_map<
             std::string,
             std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
             funcs_of_time,
         const std::true_type /*is_time_dependent*/) noexcept {
        point = the_map(point, t, funcs_of_time);
      })(std::get<Is>(maps_), mapped_point, time, functions_of_time,
         CoordinateMap_detail::is_map_time_dependent_t<Maps>{}));

  return tnsr::I<T, dim, TargetFrame>(std::move(mapped_point));
}

template <typename SourceFrame, typename TargetFrame, typename... Maps>
template <typename T, size_t... Is>
boost::optional<tnsr::I<
    T, CoordinateMap<SourceFrame, TargetFrame, Maps...>::dim, SourceFrame>>
CoordinateMap<SourceFrame, TargetFrame, Maps...>::inverse_impl(
    tnsr::I<T, dim, TargetFrame>&& target_point, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time,
    std::index_sequence<Is...> /*meta*/) const noexcept {
  boost::optional<std::array<T, dim>> mapped_point(
      make_array<T, dim>(std::move(target_point)));

  EXPAND_PACK_LEFT_TO_RIGHT(make_overloader(
      [](const auto& the_map, boost::optional<std::array<T, dim>>& point,
         const double /*t*/,
         const std::unordered_map<
             std::string,
             std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
         /*funcs_of_time*/,
         const std::false_type /*is_time_independent*/) noexcept {
        if (point) {
          if (LIKELY(not the_map.is_identity())) {
            point = the_map.inverse(point.get());
          }
        }
      },
      [](const auto& the_map, boost::optional<std::array<T, dim>>& point,
         const double t,
         const std::unordered_map<
             std::string,
             std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
             funcs_of_time,
         const std::true_type /*is_time_dependent*/) noexcept {
        if (point) {
          point = the_map.inverse(point.get(), t, funcs_of_time);
        }
        // this is the inverse function, so the iterator sequence below is
        // reversed
      })(std::get<sizeof...(Maps) - 1 - Is>(maps_), mapped_point, time,
         functions_of_time,
         CoordinateMap_detail::is_map_time_dependent_t<decltype(
             std::get<sizeof...(Maps) - 1 - Is>(maps_))>{}));

  return mapped_point
             ? tnsr::I<T, dim, SourceFrame>(std::move(mapped_point.get()))
             : boost::optional<tnsr::I<T, dim, SourceFrame>>{};
}

// define type-trait to check for time-dependent jacobian
namespace CoordinateMap_detail {
CREATE_IS_CALLABLE(jacobian)
template <typename Map, typename T>
using is_jacobian_time_dependent_t =
    CoordinateMap_detail::is_jacobian_callable_t<
        Map, std::array<std::decay_t<T>, std::decay_t<Map>::dim>, double,
        std::unordered_map<
            std::string,
            std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>>;
}  // namespace CoordinateMap_detail

template <typename SourceFrame, typename TargetFrame, typename... Maps>
template <typename T>
auto CoordinateMap<SourceFrame, TargetFrame, Maps...>::inv_jacobian_impl(
    tnsr::I<T, dim, SourceFrame>&& source_point, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const noexcept
    -> InverseJacobian<T, dim, SourceFrame, TargetFrame> {
  std::array<T, dim> mapped_point = make_array<T, dim>(std::move(source_point));

  InverseJacobian<T, dim, SourceFrame, TargetFrame> inv_jac{};

  tuple_transform(
      maps_,
      [&inv_jac, &mapped_point, time, &functions_of_time](
          const auto& map, auto index, const std::tuple<Maps...>& maps) {
        constexpr const size_t count = decltype(index)::value;

        tnsr::Ij<T, dim, Frame::NoFrame> temp_inv_jac{};

        // chooses the correct call based on time-dependence of jacobian
        auto inv_jac_overload = make_overloader(
            [](const gsl::not_null<tnsr::Ij<T, dim, Frame::NoFrame>*> t_inv_jac,
               const auto& the_map, const std::array<T, dim>& point,
               const double /*t*/,
               const std::unordered_map<
                   std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
               /*funcs_of_time*/,
               const std::false_type /*is_time_independent*/) noexcept {
              if (LIKELY(not the_map.is_identity())) {
                *t_inv_jac = the_map.inv_jacobian(point);
              } else {
                *t_inv_jac = identity<dim>(point[0]);
              }
              return nullptr;
            },
            [](const gsl::not_null<tnsr::Ij<T, dim, Frame::NoFrame>*> t_inv_jac,
               const auto& the_map, const std::array<T, dim>& point,
               const double t,
               const std::unordered_map<
                   std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
                   funcs_of_time,
               const std::true_type /*is_time_dependent*/) noexcept {
              *t_inv_jac = the_map.inv_jacobian(point, t, funcs_of_time);
              return nullptr;
            });

        if (LIKELY(count != 0)) {
          const auto& map_in_loop =
              std::get<(count != 0 ? count - 1 : 0)>(maps);
          if (LIKELY(not map_in_loop.is_identity())) {
            CoordinateMap_detail::apply_map(
                make_not_null(&mapped_point), map_in_loop, time,
                functions_of_time,
                CoordinateMap_detail::is_map_time_dependent_t<decltype(
                    map_in_loop)>{});
            inv_jac_overload(&temp_inv_jac, map, mapped_point, time,
                             functions_of_time,
                             CoordinateMap_detail::is_jacobian_time_dependent_t<
                                 decltype(map), T>{});
            std::array<T, dim> temp{};
            for (size_t source = 0; source < dim; ++source) {
              for (size_t target = 0; target < dim; ++target) {
                gsl::at(temp, target) =
                    inv_jac.get(source, 0) * temp_inv_jac.get(0, target);
                for (size_t dummy = 1; dummy < dim; ++dummy) {
                  gsl::at(temp, target) += inv_jac.get(source, dummy) *
                                           temp_inv_jac.get(dummy, target);
                }
              }
              for (size_t target = 0; target < dim; ++target) {
                inv_jac.get(source, target) = std::move(gsl::at(temp, target));
              }
            }
          }
        } else {
          inv_jac_overload(
              &temp_inv_jac, map, mapped_point, time, functions_of_time,
              CoordinateMap_detail::is_jacobian_time_dependent_t<decltype(map),
                                                                 T>{});
          for (size_t source = 0; source < dim; ++source) {
            for (size_t target = 0; target < dim; ++target) {
              inv_jac.get(source, target) =
                  std::move(temp_inv_jac.get(source, target));
            }
          }
        }
      },
      maps_);
  return inv_jac;
}

template <typename SourceFrame, typename TargetFrame, typename... Maps>
template <typename T>
auto CoordinateMap<SourceFrame, TargetFrame, Maps...>::jacobian_impl(
    tnsr::I<T, dim, SourceFrame>&& source_point, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const noexcept
    -> Jacobian<T, dim, SourceFrame, TargetFrame> {
  std::array<T, dim> mapped_point = make_array<T, dim>(std::move(source_point));
  Jacobian<T, dim, SourceFrame, TargetFrame> jac{};

  tuple_transform(
      maps_,
      [&jac, &mapped_point, time, &functions_of_time](
          const auto& map, auto index,
          const std::tuple<Maps...>& maps) noexcept {
        constexpr const size_t count = decltype(index)::value;

        tnsr::Ij<T, dim, Frame::NoFrame> noframe_jac{};

        // chooses the correct call based on time-dependence of jacobian
        auto jac_overload = make_overloader(
            [](const gsl::not_null<tnsr::Ij<T, dim, Frame::NoFrame>*>
                   no_frame_jac,
               const auto& the_map, const std::array<T, dim>& point,
               const double /*t*/,
               const std::unordered_map<
                   std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
               /*funcs_of_time*/,
               const std::false_type /*is_time_independent*/) noexcept {
              if (LIKELY(not the_map.is_identity())) {
                *no_frame_jac = the_map.jacobian(point);
              } else {
                *no_frame_jac = identity<dim>(point[0]);
              }
              return nullptr;
            },
            [](const gsl::not_null<tnsr::Ij<T, dim, Frame::NoFrame>*>
                   no_frame_jac,
               const auto& the_map, const std::array<T, dim>& point,
               const double t,
               const std::unordered_map<
                   std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
                   funcs_of_time,
               const std::true_type /*is_time_dependent*/) noexcept {
              *no_frame_jac = the_map.jacobian(point, t, funcs_of_time);
              return nullptr;
            });

        if (LIKELY(count != 0)) {
          const auto& map_in_loop =
              std::get<(count != 0 ? count - 1 : 0)>(maps);
          if (LIKELY(not map_in_loop.is_identity())) {
            CoordinateMap_detail::apply_map(
                make_not_null(&mapped_point), map_in_loop, time,
                functions_of_time,
                CoordinateMap_detail::is_map_time_dependent_t<decltype(
                    map_in_loop)>{});
            jac_overload(&noframe_jac, map, mapped_point, time,
                         functions_of_time,
                         CoordinateMap_detail::is_jacobian_time_dependent_t<
                             decltype(map), T>{});
            std::array<T, dim> temp{};
            for (size_t source = 0; source < dim; ++source) {
              for (size_t target = 0; target < dim; ++target) {
                gsl::at(temp, target) =
                    noframe_jac.get(target, 0) * jac.get(0, source);
                for (size_t dummy = 1; dummy < dim; ++dummy) {
                  gsl::at(temp, target) +=
                      noframe_jac.get(target, dummy) * jac.get(dummy, source);
                }
              }
              for (size_t target = 0; target < dim; ++target) {
                jac.get(target, source) = std::move(gsl::at(temp, target));
              }
            }
          }
        } else {
          jac_overload(
              &noframe_jac, map, mapped_point, time, functions_of_time,
              CoordinateMap_detail::is_jacobian_time_dependent_t<decltype(map),
                                                                 T>{});
          for (size_t target = 0; target < dim; ++target) {
            for (size_t source = 0; source < dim; ++source) {
              jac.get(target, source) =
                  std::move(noframe_jac.get(target, source));
            }
          }
        }
      },
      maps_);
  return jac;
}

template <typename SourceFrame, typename TargetFrame, typename... Maps>
bool operator!=(
    const CoordinateMap<SourceFrame, TargetFrame, Maps...>& lhs,
    const CoordinateMap<SourceFrame, TargetFrame, Maps...>& rhs) noexcept {
  return not(lhs == rhs);
}

template <typename SourceFrame, typename TargetFrame, typename... Maps>
auto make_coordinate_map(Maps&&... maps) noexcept
    -> CoordinateMap<SourceFrame, TargetFrame, std::decay_t<Maps>...> {
  return CoordinateMap<SourceFrame, TargetFrame, std::decay_t<Maps>...>(
      std::forward<Maps>(maps)...);
}

template <typename SourceFrame, typename TargetFrame, typename... Maps>
auto make_coordinate_map_base(Maps&&... maps) noexcept
    -> std::unique_ptr<CoordinateMapBase<
        SourceFrame, TargetFrame,
        CoordinateMap<SourceFrame, TargetFrame, std::decay_t<Maps>...>::dim>> {
  return std::make_unique<
      CoordinateMap<SourceFrame, TargetFrame, std::decay_t<Maps>...>>(
      std::forward<Maps>(maps)...);
}

template <typename SourceFrame, typename TargetFrame, typename Arg0,
          typename... Args>
auto make_vector_coordinate_map_base(Arg0&& arg_0,
                                     Args&&... remaining_args) noexcept
    -> std::vector<std::unique_ptr<
        CoordinateMapBase<SourceFrame, TargetFrame, std::decay_t<Arg0>::dim>>> {
  std::vector<std::unique_ptr<
      CoordinateMapBase<SourceFrame, TargetFrame, std::decay_t<Arg0>::dim>>>
      return_vector;
  return_vector.reserve(sizeof...(Args) + 1);
  return_vector.emplace_back(make_coordinate_map_base<SourceFrame, TargetFrame>(
      std::forward<Arg0>(arg_0)));
  EXPAND_PACK_LEFT_TO_RIGHT(return_vector.emplace_back(
      make_coordinate_map_base<SourceFrame, TargetFrame>(
          std::forward<Args>(remaining_args))));
  return return_vector;
}

template <typename SourceFrame, typename TargetFrame, size_t Dim, typename Map,
          typename... Maps>
auto make_vector_coordinate_map_base(std::vector<Map> maps,
                                     const Maps&... remaining_maps) noexcept
    -> std::vector<
        std::unique_ptr<CoordinateMapBase<SourceFrame, TargetFrame, Dim>>> {
  std::vector<std::unique_ptr<CoordinateMapBase<SourceFrame, TargetFrame, Dim>>>
      return_vector;
  return_vector.reserve(sizeof...(Maps) + 1);
  for (auto& map : maps) {
    return_vector.emplace_back(
        make_coordinate_map_base<SourceFrame, TargetFrame>(std::move(map),
                                                           remaining_maps...));
  }
  return return_vector;
}

template <typename NewMap, typename SourceFrame, typename TargetFrame,
          typename... Maps, size_t... Is>
CoordinateMap<SourceFrame, TargetFrame, Maps..., NewMap> push_back_impl(
    CoordinateMap<SourceFrame, TargetFrame, Maps...>&& old_map, NewMap new_map,
    std::index_sequence<Is...> /*meta*/) noexcept {
  return CoordinateMap<SourceFrame, TargetFrame, Maps..., NewMap>{
      std::move(std::get<Is>(old_map.maps_))..., std::move(new_map)};
}

template <typename NewMap, typename SourceFrame, typename TargetFrame,
          typename... Maps, size_t... Is>
CoordinateMap<SourceFrame, TargetFrame, NewMap, Maps...> push_front_impl(
    CoordinateMap<SourceFrame, TargetFrame, Maps...>&& old_map, NewMap new_map,
    std::index_sequence<Is...> /*meta*/) noexcept {
  return CoordinateMap<SourceFrame, TargetFrame, Maps..., NewMap>{
      std::move(new_map), std::move(std::get<Is>(old_map.maps_))...};
}

template <typename SourceFrame, typename TargetFrame, typename... Maps,
          typename NewMap>
CoordinateMap<SourceFrame, TargetFrame, Maps..., NewMap> push_back(
    CoordinateMap<SourceFrame, TargetFrame, Maps...> old_map,
    NewMap new_map) noexcept {
  return push_back_impl(std::move(old_map), std::move(new_map),
                        std::make_index_sequence<sizeof...(Maps)>{});
}

template <typename SourceFrame, typename TargetFrame, typename... Maps,
          typename NewMap>
CoordinateMap<SourceFrame, TargetFrame, NewMap, Maps...> push_front(
    CoordinateMap<SourceFrame, TargetFrame, Maps...> old_map,
    NewMap new_map) noexcept {
  return push_front_impl(std::move(old_map), std::move(new_map),
                         std::make_index_sequence<sizeof...(Maps)>{});
}
}  // namespace domain
/// \endcond
