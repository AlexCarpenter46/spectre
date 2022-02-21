# Distributed under the MIT License.
# See LICENSE.txt for details.

find_package(GSL REQUIRED)

message(STATUS "GSL libs: ${GSL_LIBRARIES}")
message(STATUS "GSL incl: ${GSL_INCLUDE_DIR}")
message(STATUS "GSL vers: ${GSL_VERSION}")

file(APPEND
  "${CMAKE_BINARY_DIR}/BuildInfo.txt"
  "GSL version: ${GSL_VERSION}\n"
  )

# Link external BLAS library. We don't need the GSL::gslcblas target.
target_link_libraries(
  GSL::gsl
  INTERFACE
  Blas
  )

set_property(
  GLOBAL APPEND PROPERTY SPECTRE_THIRD_PARTY_LIBS
  GSL::gsl GSL::gslcblas
  )
