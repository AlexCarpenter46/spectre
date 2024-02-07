// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class Informer.

#pragma once

#include <string>

/// \cond
class CkArgMsg;
/// \endcond

/// \ingroup LoggingGroup
/// The Informer manages textual output regarding the status of a simulation.
class Informer {
 public:
  /// Useful information to print at the beginning of a simulation.
  ///
  /// This includes the command used to start the executable such as
  ///
  /// ```
  /// ./MyExecutable --input-file MyInputFile.yaml
  /// ```
  ///
  /// If you used charmrun, mpirun, or something similar to start your
  /// executable, you'll only see the options that have to do with the
  /// executable itself. Meaning, for this command
  ///
  /// ```
  /// mpirun -np 4 MyExecutable --input-file MyInputFile.yaml
  /// ```
  ///
  /// only `MyExecutable` and onwards will be printed.
  static std::string startup_info(CkArgMsg* msg);

  /// Useful information to print at the end of a simulation.
  static std::string exit_info();
};
