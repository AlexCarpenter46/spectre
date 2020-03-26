# Distributed under the MIT License.
# See LICENSE.txt for details.

find_path(
    BLAZE_INCLUDE_DIR
    PATH_SUFFIXES include
    NAMES blaze/Blaze.h
    HINTS ${BLAZE_ROOT} ENV BLAZE_ROOT
    DOC "Blaze include directory. Used BLAZE_ROOT to set a search dir."
)

set(BLAZE_INCLUDE_DIRS ${BLAZE_INCLUDE_DIR})

# Extract version info from header
file(READ
  "${BLAZE_INCLUDE_DIRS}/blaze/system/Version.h"
  BLAZE_FIND_HEADER_CONTENTS)

string(REGEX MATCH "#define BLAZE_MAJOR_VERSION [0-9]+"
  BLAZE_MAJOR_VERSION "${BLAZE_FIND_HEADER_CONTENTS}")
string(REPLACE "#define BLAZE_MAJOR_VERSION " ""
  BLAZE_MAJOR_VERSION
  ${BLAZE_MAJOR_VERSION})

string(REGEX MATCH "#define BLAZE_MINOR_VERSION [0-9]+"
  BLAZE_MINOR_VERSION "${BLAZE_FIND_HEADER_CONTENTS}")
string(REPLACE "#define BLAZE_MINOR_VERSION " ""
  BLAZE_MINOR_VERSION
  ${BLAZE_MINOR_VERSION})

set(BLAZE_VERSION
  "${BLAZE_MAJOR_VERSION}.${BLAZE_MINOR_VERSION}"
  )

set(Blaze_VERSION
  "${BLAZE_MAJOR_VERSION}.${BLAZE_MINOR_VERSION}"
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Blaze
  FOUND_VAR BLAZE_FOUND
  REQUIRED_VARS BLAZE_INCLUDE_DIR BLAZE_INCLUDE_DIRS
  VERSION_VAR BLAZE_VERSION
  )
mark_as_advanced(BLAZE_INCLUDE_DIR BLAZE_INCLUDE_DIRS
  BLAZE_VERSION BLAZE_MAJOR_VERSION BLAZE_MINOR_VERSION
  )
