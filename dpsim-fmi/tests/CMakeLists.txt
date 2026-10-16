# SPDX-License-Identifier: MPL-2.0
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
	set(FMI_ARCH aarch64)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
	set(FMI_ARCH x86_64)
else()
	message(FATAL_ERROR "FMI test binaries do not support this architecture")
endif()

if(APPLE)
	set(FMI_PLATFORM "${FMI_ARCH}-darwin")
	set(FMI_SUFFIX .dylib)
elseif(WIN32)
	set(FMI_PLATFORM "${FMI_ARCH}-windows")
	set(FMI_SUFFIX .dll)
else()
	set(FMI_PLATFORM "${FMI_ARCH}-linux")
	set(FMI_SUFFIX .so)
endif()

add_executable(FMI_ImportTest ImportTest.cpp)
target_link_libraries(FMI_ImportTest PRIVATE dpsim-fmi)
target_include_directories(FMI_ImportTest PRIVATE ../examples)

set(FIXTURES linear nonlinear affine timevarying event initerror singleton missingfunction)
set(KIND 0)

foreach(FIXTURE IN LISTS FIXTURES)
	set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/${FIXTURE}")
	set(SINGLETON false)

	if(FIXTURE STREQUAL "singleton")
		set(SINGLETON true)
	endif()

	file(MAKE_DIRECTORY "${STAGE}")
	configure_file(modelDescription.xml.in "${STAGE}/modelDescription.xml" @ONLY)

	add_library(fmi_fixture_${FIXTURE} MODULE TestModel.cpp)
	target_include_directories(fmi_fixture_${FIXTURE} PRIVATE ../third_party/fmi3)
	target_compile_definitions(fmi_fixture_${FIXTURE} PRIVATE FIXTURE_KIND=${KIND})

	set_target_properties(fmi_fixture_${FIXTURE} PROPERTIES
		PREFIX ""
		OUTPUT_NAME TestModel
		SUFFIX "${FMI_SUFFIX}"
		LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
		RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	)

	set(FMU "${CMAKE_CURRENT_BINARY_DIR}/${FIXTURE}.fmu")

	add_custom_command(OUTPUT "${FMU}"
		COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip --
			modelDescription.xml binaries
		WORKING_DIRECTORY "${STAGE}"
		DEPENDS fmi_fixture_${FIXTURE} "${STAGE}/modelDescription.xml"
	)

	add_custom_target(fmi_package_${FIXTURE} DEPENDS "${FMU}")
	add_dependencies(FMI_ImportTest fmi_package_${FIXTURE})
	add_test(NAME fmi_${FIXTURE} COMMAND FMI_ImportTest "${FMU}" "${FIXTURE}")

	math(EXPR KIND "${KIND}+1")
endforeach()

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME fmi_archive_validation COMMAND ${Python3_EXECUTABLE}
	${CMAKE_CURRENT_SOURCE_DIR}/ArchiveTest.py
	$<TARGET_FILE:FMI_ImportTest>
	${CMAKE_CURRENT_BINARY_DIR}/linear.fmu
)

set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/threephase")
file(MAKE_DIRECTORY "${STAGE}")
configure_file(modelDescription3Ph.xml "${STAGE}/modelDescription.xml" COPYONLY)

add_library(fmi_fixture_threephase MODULE TestModel3Ph.cpp)
target_include_directories(fmi_fixture_threephase PRIVATE ../third_party/fmi3)

set_target_properties(fmi_fixture_threephase PROPERTIES
	PREFIX ""
	OUTPUT_NAME TestModel3Ph
	SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
)

set(FMU "${CMAKE_CURRENT_BINARY_DIR}/threephase.fmu")

add_custom_command(OUTPUT "${FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_threephase "${STAGE}/modelDescription.xml"
)

add_custom_target(fmi_package_threephase DEPENDS "${FMU}")

add_executable(FMI_ImportTest3Ph ImportTest3Ph.cpp)
target_link_libraries(FMI_ImportTest3Ph PRIVATE dpsim-fmi)
target_include_directories(FMI_ImportTest3Ph PRIVATE ../examples)
add_dependencies(FMI_ImportTest3Ph fmi_package_threephase)

add_test(NAME fmi_threephase COMMAND FMI_ImportTest3Ph "${FMU}")
add_test(NAME fmi_threephase_floating COMMAND FMI_ImportTest3Ph "${FMU}" floating)

set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/nonlinear_threephase")
file(MAKE_DIRECTORY "${STAGE}")
configure_file(modelDescription3Ph.xml "${STAGE}/modelDescription.xml" COPYONLY)

add_library(fmi_fixture_nonlinear_threephase MODULE TestModel3Ph.cpp)
target_include_directories(fmi_fixture_nonlinear_threephase PRIVATE ../third_party/fmi3)
target_compile_definitions(fmi_fixture_nonlinear_threephase PRIVATE FIXTURE_NONLINEAR=1)

set_target_properties(fmi_fixture_nonlinear_threephase PROPERTIES
	PREFIX ""
	OUTPUT_NAME TestModel3Ph
	SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
)

set(FMU "${CMAKE_CURRENT_BINARY_DIR}/nonlinear_threephase.fmu")

add_custom_command(OUTPUT "${FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_nonlinear_threephase "${STAGE}/modelDescription.xml"
)

add_custom_target(fmi_package_nonlinear_threephase DEPENDS "${FMU}")

add_executable(FMI_NonlinearTest NonlinearTest.cpp)
target_link_libraries(FMI_NonlinearTest PRIVATE dpsim-fmi)
target_include_directories(FMI_NonlinearTest PRIVATE ../examples)
add_dependencies(FMI_NonlinearTest fmi_package_nonlinear_threephase)

add_test(NAME fmi_nonlinear_threephase COMMAND FMI_NonlinearTest "${FMU}")
add_test(NAME fmi_nonlinear_reject_unconverged COMMAND FMI_NonlinearTest "${FMU}" fail)

# A completed-step event changes terminal current without reinitializing states.
# It must fail explicitly until instantaneous algebraic event iteration exists.
set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/event_current")
file(MAKE_DIRECTORY "${STAGE}")
configure_file(modelDescription3Ph.xml "${STAGE}/modelDescription.xml" COPYONLY)

add_library(fmi_fixture_event_current MODULE TestModel3Ph.cpp)
target_include_directories(fmi_fixture_event_current PRIVATE ../third_party/fmi3)
target_compile_definitions(
	fmi_fixture_event_current
	PRIVATE FIXTURE_NONLINEAR=1 FIXTURE_EVENT_JUMP=1
)

set_target_properties(fmi_fixture_event_current PROPERTIES
	PREFIX ""
	OUTPUT_NAME TestModel3Ph
	SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
)

set(FMU "${CMAKE_CURRENT_BINARY_DIR}/event_current.fmu")

add_custom_command(OUTPUT "${FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_event_current "${STAGE}/modelDescription.xml"
)

add_custom_target(fmi_package_event_current DEPENDS "${FMU}")

add_executable(FMI_EventCurrentTest EventCurrentTest.cpp)
target_link_libraries(FMI_EventCurrentTest PRIVATE dpsim-fmi)
target_include_directories(FMI_EventCurrentTest PRIVATE ../examples)
add_dependencies(FMI_EventCurrentTest fmi_package_event_current)

add_test(NAME fmi_reject_event_current_jump COMMAND FMI_EventCurrentTest "${FMU}")

# An explicitly enabled state event wraps a periodic continuous state while its
# terminal current remains continuous. The reset must become DPsim history.
set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/state_event")
file(MAKE_DIRECTORY "${STAGE}")
configure_file(modelDescriptionStateEvent.xml "${STAGE}/modelDescription.xml" COPYONLY)

add_library(fmi_fixture_state_event MODULE StateEventModel.cpp)
target_include_directories(fmi_fixture_state_event PRIVATE ../third_party/fmi3)

set_target_properties(fmi_fixture_state_event PROPERTIES
	PREFIX ""
	OUTPUT_NAME StateEventModel
	SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
)

set(FMU "${CMAKE_CURRENT_BINARY_DIR}/state_event.fmu")

add_custom_command(OUTPUT "${FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_state_event "${STAGE}/modelDescription.xml"
)

add_custom_target(fmi_package_state_event DEPENDS "${FMU}")

add_executable(FMI_StateEventTest StateEventTest.cpp)
target_link_libraries(FMI_StateEventTest PRIVATE dpsim-fmi)
target_include_directories(FMI_StateEventTest PRIVATE ../examples)
add_dependencies(FMI_StateEventTest fmi_package_state_event)

add_test(NAME fmi_state_event_reset COMMAND FMI_StateEventTest "${FMU}")

set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/sampled_threephase")
file(MAKE_DIRECTORY "${STAGE}")

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/modelDescription3Ph.xml" SAMPLED_XML)

string(
	REPLACE
		"</ModelVariables>"
		"<Float64 name=\"held_gain\" valueReference=\"22\" variability=\"discrete\" causality=\"local\" initial=\"calculated\"/><Float64 name=\"command\" valueReference=\"23\" causality=\"input\" variability=\"continuous\" start=\"1\"/></ModelVariables>"
		SAMPLED_XML
		"${SAMPLED_XML}"
)

file(WRITE "${STAGE}/modelDescription.xml" "${SAMPLED_XML}")

add_library(fmi_fixture_sampled_threephase MODULE TestModel3Ph.cpp)
target_include_directories(fmi_fixture_sampled_threephase PRIVATE ../third_party/fmi3)
target_compile_definitions(fmi_fixture_sampled_threephase PRIVATE FIXTURE_NONLINEAR=1 FIXTURE_SAMPLED=1)

set_target_properties(fmi_fixture_sampled_threephase PROPERTIES
	PREFIX ""
	OUTPUT_NAME TestModel3Ph
	SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
)

set(FMU "${CMAKE_CURRENT_BINARY_DIR}/sampled_threephase.fmu")

add_custom_command(OUTPUT "${FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_sampled_threephase "${STAGE}/modelDescription.xml"
)

add_custom_target(fmi_package_sampled_threephase DEPENDS "${FMU}")

add_executable(FMI_SampledTest SampledTest.cpp)
target_link_libraries(FMI_SampledTest PRIVATE dpsim-fmi)
target_include_directories(FMI_SampledTest PRIVATE ../examples)
add_dependencies(FMI_SampledTest fmi_package_sampled_threephase)

add_test(NAME fmi_sampled_threephase COMMAND FMI_SampledTest "${FMU}")
add_test(NAME fmi_reject_misaligned_tick COMMAND FMI_SampledTest "${FMU}" misaligned)
add_test(NAME fmi_reject_sampled_lti COMMAND FMI_ImportTest3Ph "${FMU}" reject_sampled)

# Additional runtime lifecycle regressions share the same nonlinear equations.
function(add_fmi_lifecycle_fixture name)
	cmake_parse_arguments(ARG "SNAPSHOTS;NO_COMPLETED" "" "DEFINITIONS" ${ARGN})
	set(stage "${CMAKE_CURRENT_BINARY_DIR}/stage/${name}")
	file(MAKE_DIRECTORY "${stage}")

	file(READ "${CMAKE_CURRENT_SOURCE_DIR}/modelDescription3Ph.xml" xml)

	if(ARG_NO_COMPLETED)
		string(
			REPLACE
				"needsCompletedIntegratorStep=\"true\""
				"needsCompletedIntegratorStep=\"false\""
				xml
				"${xml}"
		)
	endif()

	if(ARG_SNAPSHOTS)
		string(
			REPLACE
				"<ModelExchange "
				"<ModelExchange canGetAndSetFMUState=\"true\" "
				xml
				"${xml}"
		)

		string(
			REPLACE
				"</ModelVariables>"
				"<Float64 name=\"evaluation_count\" valueReference=\"24\" causality=\"output\" variability=\"continuous\" initial=\"calculated\"/><Float64 name=\"history\" valueReference=\"25\" causality=\"local\" variability=\"discrete\" initial=\"calculated\"/></ModelVariables>"
				xml
				"${xml}"
		)

		string(
			REPLACE
				"</ModelStructure>"
				"<Output valueReference=\"24\"/><InitialUnknown valueReference=\"24\"/><InitialUnknown valueReference=\"25\"/></ModelStructure>"
				xml
				"${xml}"
		)
	endif()

	file(WRITE "${stage}/modelDescription.xml" "${xml}")

	add_library(fmi_fixture_${name} MODULE TestModel3Ph.cpp)
	target_include_directories(fmi_fixture_${name} PRIVATE ../third_party/fmi3)
	target_compile_definitions(fmi_fixture_${name} PRIVATE FIXTURE_NONLINEAR=1 ${ARG_DEFINITIONS})

	set_target_properties(fmi_fixture_${name} PROPERTIES
		PREFIX ""
		OUTPUT_NAME TestModel3Ph
		SUFFIX "${FMI_SUFFIX}"
		LIBRARY_OUTPUT_DIRECTORY "${stage}/binaries/${FMI_PLATFORM}"
		RUNTIME_OUTPUT_DIRECTORY "${stage}/binaries/${FMI_PLATFORM}"
	)

	set(fmu "${CMAKE_CURRENT_BINARY_DIR}/${name}.fmu")

	add_custom_command(OUTPUT "${fmu}"
		COMMAND ${CMAKE_COMMAND} -E tar cf "${fmu}" --format=zip -- modelDescription.xml binaries
		WORKING_DIRECTORY "${stage}"
		DEPENDS fmi_fixture_${name} "${stage}/modelDescription.xml"
	)

	add_custom_target(fmi_package_${name} DEPENDS "${fmu}")
	add_dependencies(FMI_NonlinearTest fmi_package_${name})
	add_test(NAME fmi_${name} COMMAND FMI_NonlinearTest "${fmu}")
endfunction()

add_fmi_lifecycle_fixture(snapshot_isolation SNAPSHOTS DEFINITIONS FIXTURE_SNAPSHOTS=1)
add_fmi_lifecycle_fixture(skip_completed_step NO_COMPLETED DEFINITIONS FIXTURE_FORBID_COMPLETED=1)
