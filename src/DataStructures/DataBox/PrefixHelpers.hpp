// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
template <typename TagsList>
class Variables;
namespace db {
struct PrefixTag;
struct SimpleTag;
}  // namespace db

namespace Tags {
template <typename TagsList>
struct Variables;
}  // namespace Tags
/// \endcond

namespace db {

/// \ingroup DataBoxTagsGroup
/// \brief Create a new `tmpl::list` of tags by wrapping each tag in `TagList`
/// in `Wrapper<_, Args...>`.
template <template <typename...> class Wrapper, typename TagList,
          typename... Args>
using wrap_tags_in =
    tmpl::transform<TagList, tmpl::bind<Wrapper, tmpl::_1, tmpl::pin<Args>...>>;

namespace detail {
template <template <typename...> class Prefix, typename Tag, typename... Args>
struct add_tag_prefix_impl {
  using type = Prefix<Tag, Args...>;
};

template <template <typename...> class Prefix, typename TagList,
          typename... Args>
struct add_tag_prefix_impl<Prefix, Tags::Variables<TagList>, Args...> {
  using type = Tags::Variables<wrap_tags_in<Prefix, TagList, Args...>>;
};
}  // namespace detail

/// \ingroup DataBoxTagsGroup
/// Wrap `Tag` in `Prefix<_, Args...>`, unless `Tag` is a Tags::Variables,
/// in which case this creates a new Tags::Variables, wrapping each tag in
/// `Tag::tags_list` with `Prefix<_, Args...>`.
template <template <typename...> class Prefix, typename Tag, typename... Args>
using add_tag_prefix =
    typename detail::add_tag_prefix_impl<Prefix, Tag, Args...>::type;

namespace detail {
template <typename>
struct remove_tag_prefix_impl;

template <typename WrappedTag, template <typename...> class Prefix,
          typename... Args>
struct remove_tag_prefix_impl<Prefix<WrappedTag, Args...>> {
  using type = WrappedTag;
};

template <typename TagList>
struct remove_tag_prefix_impl<Tags::Variables<TagList>> {
  using type = Tags::Variables<
      tmpl::transform<TagList, remove_tag_prefix_impl<tmpl::_1>>>;
};
}  // namespace detail

/// \ingroup DataBoxTagsGroup
/// Remove the outer prefix from a prefixed tag `Tag`, or remove the outer
/// prefix of each tag in `Tag::tags_list` if `Tag` is a Tags::Variables.
template <typename Tag>
using remove_tag_prefix = typename detail::remove_tag_prefix_impl<Tag>::type;

namespace detail {

template <typename Tag, typename = std::nullptr_t>
struct remove_all_prefixes_impl {
  using type = Tag;
};

template <typename Tag>
struct remove_all_prefixes_impl<
    Tag, Requires<std::is_base_of_v<db::PrefixTag, Tag>>> {
  using type = typename remove_all_prefixes_impl<typename Tag::tag>::type;
};

template <typename TagList>
struct remove_all_prefixes_impl<Tags::Variables<TagList>> {
  using type = Tags::Variables<
      tmpl::transform<TagList, remove_all_prefixes_impl<tmpl::_1>>>;
};
}  // namespace detail

/// \ingroup DataBoxTagsGroup
/// Completely remove all prefix tags from a Tag, or all prefixes from
/// the tags in `Tag::tags_list` if `Tag` is a Tags::Variables.
template <typename Tag>
using remove_all_prefixes =
    typename detail::remove_all_prefixes_impl<Tag>::type;

namespace detail {
template <template <typename...> typename Wrapper, typename T, typename... Args>
struct prefix_variables {
  using type = T;
};

template <template <typename...> typename Wrapper, typename Tags,
          typename... Args>
struct prefix_variables<Wrapper, Variables<Tags>, Args...> {
  using type = Variables<::db::wrap_tags_in<Wrapper, Tags, Args...>>;
};
}  // namespace detail

/// \ingroup DataBoxTagsGroup
/// \brief Add a prefix to all tags in a Variables, leaving the
/// argument unchanged if it is not a Variables.
///
/// \see wrap_tags_in
template <template <typename...> class Wrapper, typename T, typename... Args>
using prefix_variables =
    typename detail::prefix_variables<Wrapper, T, Args...>::type;
}  // namespace db
