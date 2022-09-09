// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines which types are allowed, whether operations which certain types are
/// allowed, and other type-specific properties and configuration for
/// `TensorExpression`s
///
/// \details
/// To add support for a data type, modify the templates in this file and the
/// arithmetic operator overloads as necessary. Then, add tests as appropriate.

#pragma once

#include <complex>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Expressions/TensorExpression.hpp"
#include "DataStructures/VectorImpl.hpp"

namespace tenex {
template <typename DataType>
struct NumberAsExpression;

namespace detail {
/// \brief Whether or not `TensorExpression`s supports using a given type as a
/// numeric term
///
/// \details
/// To make it possible to use a new numeric data type as a term in
/// `TensorExpression`s, add the type to this alias and adjust other templates
/// in this file, as necessary.
///
/// \tparam X the arithmetic data type
template <typename X>
using is_supported_number_datatype =
    std::bool_constant<std::is_same_v<X, double> or
                       std::is_same_v<X, std::complex<double>>>;

/// \brief Whether or not `Tensor`s with the given data type are currently
/// supported by `TensorExpression`s
///
/// \details
/// To make it possible to use a new data type in a `Tensor` term in
/// `TensorExpression`s, add the type to this alias and adjust other templates
/// in this file, as necessary.
///
/// \tparam X the `Tensor` data type
template <typename X>
using is_supported_tensor_datatype = std::bool_constant<
    std::is_same_v<X, double> or std::is_same_v<X, std::complex<double>> or
    std::is_same_v<X, DataVector> or std::is_same_v<X, ComplexDataVector>>;

/// \brief Whether or not the given type is a `VectorImpl` type
///
/// \tparam T the given type
template <typename T>
using is_vector = std::bool_constant<std::is_base_of_v<MarkAsVectorImpl, T>>;

/// \brief If the given type is a derived `VectorImpl` type, get the base
/// `VectorImpl` type, else return the given type
///
/// \tparam T the given type
/// \tparam IsVector whether or not the given type is a `VectorImpl` type
template <typename T, bool IsVector = is_vector<T>::value>
struct upcast_if_derived_vector_type;

/// If `T` is not a `VectorImpl`, the `type` is just the input
template <typename T>
struct upcast_if_derived_vector_type<T, false> {
  using type = T;
};
/// If `T` is a `VectorImpl`, the `type` is the base `VectorImpl` type
template <typename T>
struct upcast_if_derived_vector_type<T, true> {
  using upcasted_type = typename T::BaseType;
  // if we have a derived VectorImpl, get base type, else T is a base VectorImpl
  // and we use that
  using type =
      tmpl::conditional_t<std::is_base_of_v<MarkAsVectorImpl, upcasted_type>,
                          upcasted_type, T>;
};

/// \brief Whether or not a given type is the complex-valued partner to another
/// given type
///
/// \details
/// This is used to define pairings between complex-valued types and their
/// real-valued counterparts. For example, `std::complex<double>`'s
/// real-valued partner is `double` and `ComplexDataVector`'s real-valued
/// partner is `DataVector`. Keeping track of this is useful in determining
/// which operations can and can't be performed in `TensorExpression`s.
///
/// To make `TensorExpression`s aware of a new pairing, modify a current
/// template specialization or add a new one.
///
/// \tparam MaybeComplexDataType the given type to check for being the complex
/// partner to the other type
/// \tparam OtherDataType the other type
template <typename MaybeComplexDataType, typename OtherDataType>
struct is_complex_datatype_of_impl;

/// No template specialization was matched, so it's not a known pairing
template <typename MaybeComplexDataType, typename OtherDataType>
struct is_complex_datatype_of_impl : std::false_type {};
/// std::complex<T> is the complex type to T
template <typename T>
struct is_complex_datatype_of_impl<std::complex<T>, T> : std::true_type {};
/// ComplexDataVector is the complex type to DataVector
template <>
struct is_complex_datatype_of_impl<typename ComplexDataVector::BaseType,
                                   typename DataVector::BaseType>
    : std::true_type {};

/// \brief Whether or not a given type is the complex-valued parter to another
/// given type
///
/// \details
/// See `is_complex_datatype_of_impl` for which pairings are defined
///
/// \tparam MaybeComplexDataType the given type to check for being the complex
/// partner to the other type
/// \tparam OtherDataType the other type
template <typename MaybeComplexDataType, typename OtherDataType>
using is_complex_datatype_of = is_complex_datatype_of_impl<
    typename upcast_if_derived_vector_type<MaybeComplexDataType>::type,
    typename upcast_if_derived_vector_type<OtherDataType>::type>;

/// \brief Whether or not a given type is assignable to another within
/// `TensorExpression`s
///
/// \details
/// This is used to define which types can be assigned to which when evaluating
/// the result of a `TensorExpression`. For example, you can assign a
/// `DataVector` to a `double`, but not vice versa.
///
/// To enable assignment between two types that is not yet supported, modify a
/// current template specialization or add a new one.
///
/// \tparam LhsDataType the type being assigned
/// \tparam RhsDataType the type to assign the `LhsDataType` to
template <typename LhsDataType, typename RhsDataType>
struct lhs_datatype_is_assignable_to_rhs_datatype_impl;

/// If the LHS type is the complex-valued partner to the RHS type, we can assign
/// the LHS type to the RHS type
template <typename LhsDataType, typename RhsDataType>
struct lhs_datatype_is_assignable_to_rhs_datatype_impl
    : std::bool_constant<(
          is_complex_datatype_of_impl<LhsDataType, RhsDataType>::value)> {};
/// Can assign a type to itself
template <typename X>
struct lhs_datatype_is_assignable_to_rhs_datatype_impl<X, X> : std::true_type {
};
/// Can assign a `VectorImpl` to its value type, e.g. can assign a `DataVector`
/// to a `double`
template <typename ValueType, typename VectorType>
struct lhs_datatype_is_assignable_to_rhs_datatype_impl<
    VectorImpl<ValueType, VectorType>, ValueType> : std::true_type {};
/// Can assign a complex-valued `VectorImpl` to its real component's type, e.g.
/// can assign a `ComplexDataVector` to a `double` because the underlying type
/// of `ComplexDataVector` is `std::complex<double>`, whose real component is a
/// `double`
template <typename ValueType, typename VectorType>
struct lhs_datatype_is_assignable_to_rhs_datatype_impl<
    VectorImpl<std::complex<ValueType>, VectorType>, ValueType>
    : std::true_type {};

/// \brief Whether or not a given type is assignable to another within
/// `TensorExpression`s
///
/// \details
/// See `lhs_datatype_is_assignable_to_rhs_datatype_impl` for which assignments
/// are permitted
///
/// \tparam LhsDataType the type being assigned
/// \tparam RhsDataType the type to assign the `LhsDataType` to
template <typename LhsDataType, typename RhsDataType>
struct lhs_datatype_is_assignable_to_rhs_datatype {
  using type = lhs_datatype_is_assignable_to_rhs_datatype_impl<
      typename upcast_if_derived_vector_type<LhsDataType>::type,
      typename upcast_if_derived_vector_type<RhsDataType>::type>;
};

/// \brief Get the data type of a binary operation between two data types
/// that may occur in a `TensorExpression`
///
/// \details
/// This is used to define the resulting types of binary arithmetic operations
/// within `TensorExpression`s, e.g. `double OP double = double` and
/// `double OP DataVector = DataVector`.
///
/// To enable binary operations between two types that is not yet supported,
/// modify a current template specialization or add a new one.
///
/// \tparam X1 the data type of one operand
/// \tparam X2 the data type of the other operand
template <typename X1, typename X2>
struct get_binop_datatype_impl;

/// No template specialization was matched, so it's not a known pairing
template <typename X1, typename X2>
struct get_binop_datatype_impl {
  using type = std::bool_constant<false>;
};
/// A binary operation between two terms of the same type will yield a result
/// with that type
template <typename X>
struct get_binop_datatype_impl<X, X> {
  using type = X;
};
/// A binary operation between two `VectorImpl`s of the same type will
/// yield the shared derived `VectorImpl`, e.g.
/// `DataVector OP DataVector = DataVector`
template <typename ValueType, typename VectorType>
struct get_binop_datatype_impl<VectorImpl<ValueType, VectorType>,
                               VectorImpl<ValueType, VectorType>> {
  using type = VectorType;
};
/// @{
/// A binary operation between a type `T` and `std::complex<T>` yields a
/// `std::complex<T>`
template <typename T>
struct get_binop_datatype_impl<T, std::complex<T>> {
  using type = std::complex<T>;
};
template <typename T>
struct get_binop_datatype_impl<std::complex<T>, T> {
  using type = std::complex<T>;
};
/// @}
/// @{
/// A binary operation between a `VectorImpl` and its underlying value type
/// yields the `VectorImpl`, e.g. `DataVector OP double = DataVector`
template <typename ValueType, typename VectorType>
struct get_binop_datatype_impl<VectorImpl<ValueType, VectorType>, ValueType> {
  using type = VectorType;
};
template <typename ValueType, typename VectorType>
struct get_binop_datatype_impl<ValueType, VectorImpl<ValueType, VectorType>> {
  using type = VectorType;
};
/// @}
/// @{
/// A binary operation between a complex-valued `VectorImpl` and its real
/// component's type yields the `VectorImpl`, e.g.
/// `ComplexDataVector OP double = ComplexDataVector`
template <typename ValueType, typename VectorType>
struct get_binop_datatype_impl<VectorImpl<std::complex<ValueType>, VectorType>,
                               ValueType> {
  using type = VectorType;
};
template <typename ValueType, typename VectorType>
struct get_binop_datatype_impl<
    ValueType, VectorImpl<std::complex<ValueType>, VectorType>> {
  using type = VectorType;
};
/// @}
// HACKS
template <typename X1>
struct get_complex_datatype;
template <typename X1>
struct get_complex_datatype {
  using type = std::bool_constant<false>;
};
template <>
struct get_complex_datatype<double> {
  using type = std::complex<double>;
};
template <>
struct get_complex_datatype<DataVector> {
  using type = ComplexDataVector;
};
// e.g. Complex DataVector = DataVector + std::complex<double>
template <typename ValueType, typename VectorType>
struct get_binop_datatype_impl<VectorImpl<ValueType, VectorType>,
                               std::complex<ValueType>> {
  using type = typename get_complex_datatype<VectorType>::type;
};
template <typename ValueType, typename VectorType>
struct get_binop_datatype_impl<std::complex<ValueType>,
                               VectorImpl<ValueType, VectorType>> {
  using type = typename get_complex_datatype<VectorType>::type;
};
/// @{
/// A binary operation between a `DataVector` and a `ComplexDataVector` yields a
/// `ComplexDataVector`
template <>
struct get_binop_datatype_impl<typename ComplexDataVector::BaseType,
                               typename DataVector::BaseType> {
  using type = ComplexDataVector;
};
template <>
struct get_binop_datatype_impl<typename DataVector::BaseType,
                               typename ComplexDataVector::BaseType> {
  using type = ComplexDataVector;
};
/// @}

/// \brief Get the data type of a binary operation between two data types
/// that may occur in a `TensorExpression`
///
/// \details
/// See `get_binop_datatype_impl` for which data type combinations have a
/// defined result type
///
/// \tparam X1 the data type of one operand
/// \tparam X2 the data type of the other operand
template <typename X1, typename X2>
struct get_binop_datatype {
  using type = typename get_binop_datatype_impl<
      typename upcast_if_derived_vector_type<X1>::type,
      typename upcast_if_derived_vector_type<X2>::type>::type;

  static_assert(
      not std::is_same_v<type, std::bool_constant<false>>,
      "You are attempting to perform a binary arithmetic operation between "
      "two data types, but the data type of the result is not known within "
      "TensorExpressions. See tenex::detail::get_binop_datatype_impl.");
};

/// \brief Whether or not it is permitted to perform binary arithmetic
/// operations with the given types within `TensorExpression`
///
/// \details
/// See `get_binop_datatype_impl` for which data type combinations have a
/// defined result type
///
/// \tparam X1 the data type of one operand
/// \tparam X2 the data type of the other operand
template <typename X1, typename X2>
struct binop_datatypes_are_supported {
  using type = std::bool_constant<not std::is_same_v<
      typename get_binop_datatype_impl<
          typename upcast_if_derived_vector_type<X1>::type,
          typename upcast_if_derived_vector_type<X2>::type>::type,
      std::bool_constant<false>>>;
};

/// \brief Whether or not it is permitted to perform binary arithmetic
/// operations with `Tensor`s with the given types within `TensorExpression`s
///
/// \details
/// This is used to define which data types can be contained by the two
/// `Tensor`s in a binary operation, e.g.
/// `Tensor<ComplexDataVector>() OP Tensor<DataVector>()` is permitted, but
/// `Tensor<DataVector>() OP Tensor<double>()` is not.
///
/// To enable binary operations between `Tensor`s with types that are not yet
/// supported, modify a current template specialization or add a new one.
///
/// \tparam X1 the data type of one `Tensor` operand
/// \tparam X2 the data type of the other `Tensor` operand
template <typename X1, typename X2>
struct tensor_binop_datatypes_are_supported_impl;

/// Can only do `Tensor<X1>() OP Tensor<X2>()` if `X1 == X2` or if `X1` and
/// `X2` are real/complex partners like `DataVector` and `ComplexDataVector`
/// (see `is_complex_datatype_of_impl`)
template <typename X1, typename X2>
struct tensor_binop_datatypes_are_supported_impl
    : std::bool_constant<(std::is_same_v<X1, X2> or
                          is_complex_datatype_of<X1, X2>::value or
                          is_complex_datatype_of<X2, X1>::value)> {};

/// \brief Whether or not it is permitted to perform binary arithmetic
/// operations with `Tensor`s with the given types within `TensorExpression`s
///
/// \details
/// See `tensor_binop_datatypes_are_supported_impl` for which data type
/// combinations are permitted
///
/// \tparam X1 the data type of one `Tensor` operand
/// \tparam X2 the data type of the other `Tensor` operand
template <typename X1, typename X2>
struct tensor_binop_datatypes_are_supported {
  static_assert(
      is_supported_tensor_datatype<X1>::value and
          is_supported_tensor_datatype<X2>::value,
      "Cannot perform binary operations between the two Tensors with the "
      "given data types because at least one of the data types is not "
      "supported by TensorExpressions. See "
      "tenex::detail::is_supported_tensor_datatype.");
  using type = typename tensor_binop_datatypes_are_supported_impl<X1, X2>::type;
};

/// \brief Whether or not it is permitted to perform binary arithmetic
/// operations with `TensorExpression`s, based on their data types
///
/// \details
/// This is used to define which data types can be contained by the two
/// `TensorExpression`s in a binary operation, e.g.
/// `Tensor<DataVector>() OP double` and
/// `Tensor<ComplexDataVector>() OP Tensor<DataVector>()` are permitted, but
/// `Tensor<DataVector>() OP Tensor<double>()` is not. This differs from
/// `tensor_binop_datatypes_are_supported` in that
/// `tensorexpression_binop_datatypes_are_supported_impl` handles all derived
/// `TensorExpression` types, whether they represent `Tensor`s or numbers.
/// `tensor_binop_datatypes_are_supported` only handles the cases where both
/// `TensorExpression`s represent `Tensor`s.
///
/// To enable binary operations between `TensorExpression`s with types that
/// are not yet supported, modify a current template specialization or add a
/// new one.
///
/// \tparam T1 the first `TensorExpression` operand
/// \tparam T2 the second `TensorExpression` operand
template <typename T1, typename T2>
struct tensorexpression_binop_datatypes_are_supported_impl;

/// Since `T1` and `T2` represent `Tensor`s, check if we can do
/// `Tensor<T1::type>() OP Tensor<T2::type>()`
template <typename T1, typename T2>
struct tensorexpression_binop_datatypes_are_supported_impl {
  using type = typename tensor_binop_datatypes_are_supported_impl<
      typename T1::type, typename T2::type>::type;
};
/// @{
/// Can do `Tensor<X>() OP NUMBER` if we can do `X OP NUMBER`
template <typename TensorExpressionType, typename NumberType>
struct tensorexpression_binop_datatypes_are_supported_impl<
    TensorExpressionType, NumberAsExpression<NumberType>> {
  static_assert(
      is_supported_tensor_datatype<
          typename TensorExpressionType::type>::value and
          is_supported_number_datatype<NumberType>::value,
      "Cannot perform binary operations between Tensor and number with the "
      "given data types because at least one of the data types is not "
      "supported by TensorExpressions. See "
      "tenex::detail::is_supported_number_datatype and "
      "tenex::detail::is_supported_tensor_datatype.");
  using type = std::bool_constant<binop_datatypes_are_supported<
      typename TensorExpressionType::type, NumberType>::type::value>;
};
template <typename NumberType, typename TensorExpressionType>
struct tensorexpression_binop_datatypes_are_supported_impl<
    NumberAsExpression<NumberType>, TensorExpressionType> {
  using type = typename tensorexpression_binop_datatypes_are_supported_impl<
      TensorExpressionType, NumberAsExpression<NumberType>>::type;
};
/// @}

/// \brief Whether or not it is permitted to perform binary arithmetic
/// operations with `TensorExpression`s, based on their data types
///
/// \details
/// See `tensorexpression_binop_datatypes_are_supported_impl` for which data
/// type combinations are permitted
///
/// \tparam T1 the first `TensorExpression` operand
/// \tparam T2 the second `TensorExpression` operand
template <typename T1, typename T2>
struct tensorexpression_binop_datatypes_are_supported {
  static_assert(
      std::is_base_of_v<Expression, T1> and std::is_base_of_v<Expression, T2>,
      "Template arguments to "
      "tenex::detail::tensorexpression_binop_datatypes_are_supported must be "
      "TensorExpressions.");
  using type =
      typename tensorexpression_binop_datatypes_are_supported_impl<T1,
                                                                   T2>::type;
};

/// \brief The maximum number of arithmetic tensor operations allowed in a
/// `TensorExpression` subtree before having it be a splitting point in the
/// overall RHS expression, according to the data type held by the `Tensor`s in
/// the expression
///
/// \details
/// To enable splitting for `TensorExpression`s with data type, define a
/// template specialization below for your data type and set the `value`.
///
/// Before defining a max operations cap for some data type, the change should
/// first be justified by benchmarking many different tensor expressions before
/// and after introducing the new cap. The optimal cap will likely be
/// hardware-dependent, so fine-tuning this would ideally involve benchmarking
/// on each hardware architecture and then controling the value based on the
/// hardware.
template <typename DataType>
struct max_num_ops_in_sub_expression_impl {
  // effectively, no splitting for any unspecialized template type
  static constexpr size_t value = std::numeric_limits<size_t>::max();
};

/// \brief When the data type of the result of a `TensorExpression` is
/// `DataVector`, the maximum number of arithmetic tensor operations allowed in
/// a subtree before having it be a splitting point in the overall RHS
/// expression
///
/// \details
/// The current value set for when the data type is `DataVector` was benchmarked
/// by compiling with clang-10 Release and running on Intel(R) Xeon(R)
/// CPU E5-2630 v4 @ 2.20GHz.
template <>
struct max_num_ops_in_sub_expression_impl<DataVector> {
  static constexpr size_t value = 8;
};

/// \brief When the data type of the result of a `TensorExpression` is
/// `ComplexDataVector`, the maximum number of arithmetic tensor operations
/// allowed in a subtree before having it be a splitting point in the overall
/// RHS expression
///
/// \details
/// The current value set for when the data type is `ComplexDataVector` is set
/// to the value for `DataVector`, but the best `value` for `ComplexDataVector`
/// should also be investigated and fine-tuned.
template <>
struct max_num_ops_in_sub_expression_impl<ComplexDataVector> {
  static constexpr size_t value =
      max_num_ops_in_sub_expression_impl<DataVector>::value;
};

/// \brief Get maximum number of arithmetic tensor operations allowed in a
/// `TensorExpression` subtree before having it be a splitting point in the
/// overall RHS expression, according to the `DataType` held by the `Tensor`s in
/// the expression
template <typename DataType>
inline constexpr size_t max_num_ops_in_sub_expression =
    max_num_ops_in_sub_expression_impl<DataType>::value;
}  // namespace detail
}  // namespace tenex
