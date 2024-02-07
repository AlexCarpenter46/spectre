// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Informer/Informer.hpp"

#include <charm++.h>
#include <charm.h>
#include <sstream>
#include <string>

#include "Informer/InfoFromBuild.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/System/ParallelInfo.hpp"

std::string Informer::startup_info(CkArgMsg* msg) {
  std::stringstream ss{};
  ss << "\n"
     << "Executing '" << executable_name() << "' using "
     << sys::number_of_procs() << " processors.\n"
     << "Launch command line: ";
  for (int i = 0; i < msg->argc - 1; i++) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    ss << msg->argv[i] << " ";
  }
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  ss << msg->argv[msg->argc - 1];

  ss << "\n"
     << "Charm++ startup time in seconds: " << sys::wall_time() << "\n"
     << "Date and time at startup: " << current_date_and_time() << "\n"
     << info_from_build() << "\n";
  return ss.str();
}

std::string Informer::exit_info() {
  std::stringstream ss{};
  ss << "\n"
     << "Done!\n"
     << "Wall time: "  << sys::pretty_wall_time() << "\n"
     << "Date and time at completion: " << current_date_and_time() << "\n";
  return ss.str();
}
