#
# CMake support layer for SuiteSparse
#
# Copyright 2016-2021 Sergiu Deitsch <sergiu.deitsch@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

cmake_minimum_required (VERSION 3.9)

include (CheckLanguage)

project (SuiteSparse LANGUAGES C CXX VERSION 5.10.1
  DESCRIPTION "A suite of sparse matrix packages"
  HOMEPAGE_URL http://faculty.cse.tamu.edu/davis/suitesparse.html
)

set (CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (CMAKE_DEBUG_POSTFIX _debug)
#set (CMAKE_C_FLAGS "-Wall -pedantic")
# set(CMAKE_C_CLANG_TIDY 
#   clang-tidy;
#   -header-filter=.;
#   -checks=-*,clang-analyzer-*,bugprone-*,google-*,performance-*;)

include (CheckCCompilerFlag)
include (CheckFunctionExists)
include (CMakeDependentOption)
include (FeatureSummary)
include (GenerateExportHeader)
include (GNUInstallDirs)

set (CMAKE_C_VISIBILITY_PRESET hidden)
set (CMAKE_CXX_VISIBILITY_PRESET hidden)
set (CMAKE_VISIBILITY_INLINES_HIDDEN 1)

check_c_compiler_flag (-fPIC HAVE_FPIC)
check_c_compiler_flag (-Qstd=c99 HAVE_QSTD_C99)

option (WITH_LGPL "Enable GNU LGPL modules" ON)
cmake_dependent_option (WITH_GPL "Enable GNU GPL modules" ON "WITH_LGPL" OFF)

# SuiteSparse version
set (MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set (MINOR_VERSION ${PROJECT_VERSION_MINOR})
set (VERSION_PATCH ${PROJECT_VERSION_PATCH})

set (SHORT_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}")

if (NOT VERSION_PATCH EQUAL 0)
  set (SHORT_VERSION "${SHORT_VERSION}.${VERSION_PATCH}")
endif (NOT VERSION_PATCH EQUAL 0)

set (VERSION ${PROJECT_VERSION})

set (CPACK_PACKAGE_NAME SuiteSparse)
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SuiteSparse: an extended version of CSparse")
set (CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION})
set (CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION})
set (CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set (CPACK_PACKAGE_VERSION ${SHORT_VERSION})

include (CMakePackageConfigHelpers)
include (CMakePushCheckState)
include (CPack)
include (CheckIncludeFile)
include (CheckLibraryExists)
include (CheckSymbolExists)
include (CheckTypeSize)
include (ComponentVersion)
include (GenerateSourceFiles)

cmake_push_check_state (RESET)

check_include_file (complex.h HAVE_COMPLEX_H)
check_include_file (fenv.h HAVE_FENV_H)
check_library_exists (m ceil "" HAVE_M)

cmake_push_check_state (RESET)
set (CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)

if (HAVE_M)
  set (CMAKE_REQUIRED_LIBRARIES m)
endif (HAVE_M)

check_symbol_exists (feenableexcept "fenv.h" HAVE_FEENABLEEXCEPT)
cmake_pop_check_state ()

check_symbol_exists (_controlfp "float.h" HAVE__CONTROLFP)

set (CMAKE_EXTRA_INCLUDE_FILES complex.h)
check_type_size (_C_double_complex HAVE_C_DOUBLE_COMPLEX BUILTIN_TYPES_ONLY)

source_group ("Generated Files" REGULAR_EXPRESSION ".+_(di|dl|zi|zl)(_.+)?\\.c$")

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_property (DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS
    _CRT_SECURE_NO_WARNINGS)
endif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

set (AMD_SRCS
  AMD/Source/amd_1.c
  AMD/Source/amd_2.c
  AMD/Source/amd_aat.c
  AMD/Source/amd_control.c
  AMD/Source/amd_defaults.c
  AMD/Source/amd_dump.c
  AMD/Source/amd_global.c
  AMD/Source/amd_info.c
  AMD/Source/amd_order.c
  AMD/Source/amd_post_tree.c
  AMD/Source/amd_postorder.c
  AMD/Source/amd_preprocess.c
  AMD/Source/amd_valid.c
)

set (AMD_HDRS
  AMD/Include/amd.h
)

set (COLAMD_SRCS
  COLAMD/Source/colamd.c
)

set (COLAMD_HDRS
  COLAMD/Include/colamd.h
)

set (BTF_SRCS
  BTF/Source/btf_maxtrans.c
  BTF/Source/btf_order.c
  BTF/Source/btf_strongcomp.c
)

set (BTF_HDRS
  BTF/Include/btf.h
)

set (KLU_SRCS
  KLU/Source/klu.c
  KLU/Source/klu_analyze.c
  KLU/Source/klu_analyze_given.c
  KLU/Source/klu_compute_path.c
  KLU/Source/klu_defaults.c
  KLU/Source/klu_diagnostics.c
  KLU/Source/klu_dump.c
  KLU/Source/klu_extract.c
  KLU/Source/klu_factor.c
  KLU/Source/klu_free_numeric.c
  KLU/Source/klu_free_symbolic.c
  KLU/Source/klu_kernel.c
  KLU/Source/klu_memory.c
  KLU/Source/klu_print.c
  KLU/Source/klu_partial_factorization_path.c
  KLU/Source/klu_partial_refactorization_restart.c
  KLU/Source/klu_refactor.c
  KLU/Source/klu_scale.c
  KLU/Source/klu_solve.c
  KLU/Source/klu_sort.c
  KLU/Source/klu_tsolve.c
)

set (KLU
  KLU/Include/klu.h
  KLU/Include/klu_version.h
)

generate_source_files (GEN_AMD_SRCS_DINT FILES "${AMD_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_AMD_SRCS_DLONG FILES "${AMD_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_COLAMD_SRCS_DINT FILES "${COLAMD_SRCS}" SUFFIX _di
COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_COLAMD_SRCS_DLONG FILES "${COLAMD_SRCS}" SUFFIX _dl
COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_BTF_SRCS_DINT FILES "${BTF_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_BTF_SRCS_DLONG FILES "${BTF_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_KLU_SRCS_DINT FILES "${KLU_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_KLU_SRCS_DLONG FILES "${KLU_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

if (HAVE_COMPLEX_H)
  #list (APPEND PROJECT_SRCS Source/cs_convert.c)
endif (HAVE_COMPLEX_H)

if (MSVC)
  add_definitions (-D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

# output directory for binaries
set (OUTPUT_DIR ${CMAKE_BINARY_DIR})

# put binaries in a single directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIR})

add_library (amd
  AMD/Include/amd_internal.h
  ${AMD_HDRS}
  ${GEN_AMD_SRCS_DINT}
  ${GEN_AMD_SRCS_DLONG}
)

add_library (colamd
  ${COLAMD_HDRS}
  ${GEN_COLAMD_SRCS_DINT}
  ${GEN_COLAMD_SRCS_DLONG}
)

add_library (btf
  BTF/Include/btf_internal.h
  ${BTF_HDRS}
  ${GEN_BTF_SRCS_DINT}
  ${GEN_BTF_SRCS_DLONG}
)

set (SUITESPARSECONFIG_HDRS
  SuiteSparse_config/SuiteSparse_config.h
)

add_library (suitesparseconfig
  SuiteSparse_config/SuiteSparse_config.c
  ${SUITESPARSECONFIG_HDRS}
)

set (KLU_HDRS
  KLU/Include/klu.h
  KLU/Include/klu_version.h
)

add_library (klu KLU/Include/klu_internal.h
  ${GEN_KLU_SRCS_DINT}
  ${GEN_KLU_SRCS_DLONG}
  ${KLU_HDRS}
)

target_link_libraries (klu PUBLIC amd btf colamd)


set (SuiteSparse_AMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/AMD)
set (SuiteSparse_BTF_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/BTF)
set (SuiteSparse_COLAMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/COLAMD)
set (SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse)
set (SuiteSparse_KLU_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/KLU)

if (HAVE_M)
  target_link_libraries (suitesparseconfig INTERFACE m)
endif (HAVE_M)

target_include_directories (suitesparseconfig PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>"
)

target_include_directories (colamd PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_COLAMD_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_COLAMD_INCLUDE_DIR}>")
target_include_directories (amd PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_AMD_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_AMD_INCLUDE_DIR}>")
target_include_directories (btf PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_BTF_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_BTF_INCLUDE_DIR}>")
target_include_directories (klu PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_KLU_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KLU/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_KLU_INCLUDE_DIR}>")

target_link_libraries (amd PUBLIC suitesparseconfig)
target_link_libraries (colamd PUBLIC suitesparseconfig)
target_link_libraries (btf PUBLIC suitesparseconfig)
target_link_libraries (klu PUBLIC suitesparseconfig)

if (NOT HAVE_COMPLEX_H)
  #target_compile_definitions (SuiteSparse PUBLIC NCOMPLEX)
endif (NOT HAVE_COMPLEX_H)

if (HAVE_C_DOUBLE_COMPLEX)
  #target_compile_definitions (SuiteSparse PUBLIC HAVE_C_DOUBLE_COMPLEX)
endif (HAVE_C_DOUBLE_COMPLEX)

set (_SuiteSparse_TARGETS amd btf klu colamd suitesparseconfig)

set_property (TARGET suitesparseconfig PROPERTY NAME_ALIAS SuiteSparse)
set_property (TARGET suitesparseconfig PROPERTY HEADERS_VARIABLE
  SUITESPARSECONFIG_HDRS)
set_property (TARGET suitesparseconfig PROPERTY INSTALL_DESTINATION
  ${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR})

foreach (_target ${_SuiteSparse_TARGETS})
  get_property (_alias TARGET ${_target} PROPERTY NAME_ALIAS)
  string (TOUPPER "${_target}" _upper_target)

  if ("${_alias}" STREQUAL "")
    set (_alias ${_upper_target})
    set (_export_prefix ${_target})
  else ("${_alias}" STREQUAL "")
    set (_export_prefix ${_alias})
  endif ("${_alias}" STREQUAL "")

  get_property (_target_include_dir TARGET ${_target} PROPERTY
    INSTALL_DESTINATION)

  if ("${_target_include_dir}" STREQUAL "")
    set (_target_include_dir ${SuiteSparse_${_alias}_INCLUDE_DIR})
  endif ("${_target_include_dir}" STREQUAL "")

  set (_export_header
    ${CMAKE_CURRENT_BINARY_DIR}/${_target_include_dir}/${_export_prefix}_export.h)

  # Generate an export header with the name <target>_export.h for each mandatory
  # target.
  generate_export_header (${_target} EXPORT_FILE_NAME ${_export_header})

  # Add the export header to the target
  target_sources (${_target} PRIVATE ${_export_header})

  get_property (_HDRS_VAR TARGET ${_target} PROPERTY HEADERS_VARIABLE)

  if ("${_HDRS_VAR}" STREQUAL "")
    set (_HDRS_VAR ${_upper_target}_HDRS)
  endif ("${_HDRS_VAR}" STREQUAL "")

  # Install headers
  install (FILES ${${_HDRS_VAR}} ${_export_header} DESTINATION
    ${_target_include_dir} COMPONENT Development)
endforeach (_target)

set (SuiteSparse_CMake_LIBDIR ${CMAKE_INSTALL_LIBDIR}/SuiteSparse)
# suitesparse-config.cmake relative install directory
set (SuiteSparse_CMake_INSTALLDIR ${SuiteSparse_CMake_LIBDIR}/cmake)

target_sources (suitesparseconfig PRIVATE ${_output_CMake_MODULES})

install (TARGETS ${_SuiteSparse_TARGETS}
  EXPORT suitesparse-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  LIBRARY DESTINATION ${SuiteSparse_CMake_LIBDIR} COMPONENT Development
  ARCHIVE DESTINATION ${SuiteSparse_CMake_LIBDIR} COMPONENT Development)

#set (SuiteSparse_NAMESPACE SuiteSparse::)

configure_package_config_file (suitesparse-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config.cmake
  INSTALL_DESTINATION ${SuiteSparse_CMake_INSTALLDIR})
write_basic_package_version_file (suitesparse-config-version.cmake
  VERSION ${VERSION} COMPATIBILITY AnyNewerVersion)

export (TARGETS ${_SuiteSparse_TARGETS} NAMESPACE ${SuiteSparse_NAMESPACE} FILE
  suitesparse-targets.cmake)
export (PACKAGE SuiteSparse)

install (EXPORT suitesparse-targets NAMESPACE ${SuiteSparse_NAMESPACE}
  DESTINATION ${SuiteSparse_CMake_INSTALLDIR})
install (FILES
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config-version.cmake
  DESTINATION ${SuiteSparse_CMake_INSTALLDIR}
  COMPONENT Development)

set_component_version (amd 2.4.6)
set_component_version (btf 1.2.6)
set_component_version (colamd 2.9.6)
set_component_version (klu 1.3.9)
set_component_version (suitesparseconfig "${VERSION}")

if (HAVE_QSTD_C99)
  target_compile_options (suitesparseconfig PUBLIC $<$<C_COMPILER_ID:Intel>:-Qstd=c99>)
endif (HAVE_QSTD_C99)

set_target_properties (amd PROPERTIES EXPORT_NAME AMD)
set_target_properties (colamd PROPERTIES EXPORT_NAME COLAMD)
set_target_properties (btf PROPERTIES EXPORT_NAME BTF)
set_target_properties (klu PROPERTIES EXPORT_NAME KLU)
set_target_properties (suitesparseconfig PROPERTIES EXPORT_NAME Config)

set_target_properties (${_SuiteSparse_TARGETS} PROPERTIES
  POSITION_INDEPENDENT_CODE ON)

feature_summary (WHAT ALL)

add_executable(klu_test_partial_factorization_path KLU/Demo/klu_test_partial_factorization_path.c)
target_link_libraries(klu_test_partial_factorization_path PRIVATE klu)
add_executable(klu_test_partial_refactorization_restart KLU/Demo/klu_test_partial_refactorization_restart.c)
target_link_libraries(klu_test_partial_refactorization_restart PRIVATE klu)

enable_testing()

add_test(
  NAME klu_test_partial_factorization_path
  COMMAND $<TARGET_FILE:klu_test_partial_factorization_path>
)
add_test(
  NAME klu_test_partial_refactorization_restart
  COMMAND $<TARGET_FILE:klu_test_partial_refactorization_restart>
)
