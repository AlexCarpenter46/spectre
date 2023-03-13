# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "PyGeneralRelativity")

spectre_python_add_module(
    GeneralRelativity
    LIBRARY_NAME ${LIBRARY}
    SOURCES
    Bindings.cpp
    Psi4Real.cpp
    PYTHON_FILES
    __init__.py
)

spectre_python_headers(
    ${LIBRARY}
    INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
    HEADERS
    Psi4Real.hpp
)

spectre_python_link_libraries(
    ${LIBRARY}
    PRIVATE
    DataStructures
    GeneralRelativity
    pybind11::module
    Utilities
)

spectre_python_add_dependencies(
  ${LIBRARY}
  PyDataStructures
  PySpectral
  PyTensor
  )
