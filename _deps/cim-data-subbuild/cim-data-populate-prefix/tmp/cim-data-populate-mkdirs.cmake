# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/dpsim/_deps/cim-data-src")
  file(MAKE_DIRECTORY "/dpsim/_deps/cim-data-src")
endif()
file(MAKE_DIRECTORY
  "/dpsim/_deps/cim-data-build"
  "/dpsim/_deps/cim-data-subbuild/cim-data-populate-prefix"
  "/dpsim/_deps/cim-data-subbuild/cim-data-populate-prefix/tmp"
  "/dpsim/_deps/cim-data-subbuild/cim-data-populate-prefix/src/cim-data-populate-stamp"
  "/dpsim/_deps/cim-data-subbuild/cim-data-populate-prefix/src"
  "/dpsim/_deps/cim-data-subbuild/cim-data-populate-prefix/src/cim-data-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/dpsim/_deps/cim-data-subbuild/cim-data-populate-prefix/src/cim-data-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/dpsim/_deps/cim-data-subbuild/cim-data-populate-prefix/src/cim-data-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
