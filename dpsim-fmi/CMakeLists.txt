# SPDX-License-Identifier: MPL-2.0
enable_language(C)
find_package(LibXml2 REQUIRED)
include(FetchContent)

FetchContent_Declare(dpsim_miniz
	GIT_REPOSITORY https://github.com/richgel999/miniz.git
	GIT_TAG 293d4db1b7d0ffee9756d035b9ac6f7431ef8492 # 3.0.2
)

# Build only the ZIP library; do not inherit miniz's project-wide options.
FetchContent_GetProperties(dpsim_miniz)
if(NOT dpsim_miniz_POPULATED)
	cmake_policy(PUSH)
	if(POLICY CMP0169)
		cmake_policy(SET CMP0169 OLD)
	endif()
	FetchContent_Populate(dpsim_miniz)
	cmake_policy(POP)
endif()

configure_file(src/miniz_export.h.in "${CMAKE_CURRENT_BINARY_DIR}/miniz_export.h" COPYONLY)

add_library(dpsim-fmi-zip STATIC
	${dpsim_miniz_SOURCE_DIR}/miniz.c
	${dpsim_miniz_SOURCE_DIR}/miniz_zip.c
	${dpsim_miniz_SOURCE_DIR}/miniz_tinfl.c
	${dpsim_miniz_SOURCE_DIR}/miniz_tdef.c
)

set_target_properties(dpsim-fmi-zip PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(dpsim-fmi-zip PRIVATE
	${dpsim_miniz_SOURCE_DIR}
	${CMAKE_CURRENT_BINARY_DIR}
)

add_library(dpsim-fmi
	src/ModelExchange.cpp
	src/LinearModel.cpp
	src/LinearSSN.cpp
	src/NonlinearSSN.cpp
)

target_include_directories(dpsim-fmi
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	PRIVATE
		third_party/fmi3
		${dpsim_miniz_SOURCE_DIR}
		${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(dpsim-fmi
	PUBLIC
		dpsim
	PRIVATE
		LibXml2::LibXml2
		dpsim-fmi-zip
		${CMAKE_DL_LIBS}
)

target_compile_options(dpsim-fmi PRIVATE ${DPSIM_CXX_FLAGS})

install(TARGETS dpsim-fmi dpsim-fmi-zip EXPORT dpsim
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/dpsim-fmi
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES third_party/fmi3/LICENSE.txt
	DESTINATION ${CMAKE_INSTALL_DATADIR}/dpsim-fmi/licenses
	RENAME FMI-LICENSE.txt
)

install(FILES ${dpsim_miniz_SOURCE_DIR}/LICENSE
	DESTINATION ${CMAKE_INSTALL_DATADIR}/dpsim-fmi/licenses
	RENAME miniz-LICENSE.txt
)

if(DPSIM_BUILD_EXAMPLES)
	add_executable(FMI_DiscreteProbe examples/FMI_DiscreteProbe.cpp)
	target_link_libraries(FMI_DiscreteProbe PRIVATE dpsim-fmi)

	add_executable(FMI_Inspect examples/FMI_Inspect.cpp)
	target_link_libraries(FMI_Inspect PRIVATE dpsim-fmi)

	add_executable(FMI_EMT_Ph3_Prescribed examples/FMI_EMT_Ph3_Prescribed.cpp)
	target_link_libraries(FMI_EMT_Ph3_Prescribed PRIVATE dpsim-fmi)

	add_executable(FMI_EMT_Ph3_Nonlinear examples/FMI_EMT_Ph3_Nonlinear.cpp)
	target_link_libraries(FMI_EMT_Ph3_Nonlinear PRIVATE dpsim-fmi)

	add_executable(FMI_EMT_Ph3_LinearSSN examples/FMI_EMT_Ph3_LinearSSN.cpp)
	target_link_libraries(FMI_EMT_Ph3_LinearSSN PRIVATE dpsim-fmi)

	add_executable(FMI_EMT_LinearSSN examples/FMI_EMT_LinearSSN.cpp)
	target_link_libraries(FMI_EMT_LinearSSN PRIVATE dpsim-fmi)
endif()

option(DPSIM_FMI_BUILD_TESTS "Build self-contained FMI importer tests" OFF)

if(DPSIM_FMI_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
