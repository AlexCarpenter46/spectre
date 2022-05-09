// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class h5::Object abstract base class

#pragma once

namespace h5 {
/*!
 * \ingroup HDF5Group
 * \brief Abstract base class representing an object in an HDF5 file
 */
class Object {
 public:
  /// \cond HIDDEN_SYMBOLS
  Object() = default;
  Object(const Object& /*rhs*/) = delete;
  Object& operator=(const Object& /*rhs*/) = delete;
  Object(Object&& /*rhs*/) = delete;             // NOLINT
  Object& operator=(Object&& /*rhs*/) = delete;  // NOLINT
  virtual ~Object() = default;
  /// \endcond

  /// Return the path to the subfile where this object is stored
  virtual const std::string& subfile_path() const = 0;
};
}  // namespace h5
