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
target_include_directories(FMI_ImportTest PRIVATE ../examples/common)

set(FIXTURES linear nonlinear affine timevarying event initerror singleton missingfunction)
set(KIND 0)

foreach(FIXTURE IN LISTS FIXTURES)
	set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/${FIXTURE}")
	set(SINGLETON false)
	set(EVENT_INDICATOR "")
	if(FIXTURE STREQUAL "event")
		set(EVENT_INDICATOR "<EventIndicator valueReference=\"4\"/>")
	endif()

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

# FMI 2.0 Model Exchange compatibility fixture.  It intentionally exercises
# the same scalar linear importer path as the FMI 3.0 fixture.
set(FMI2_PLATFORM linux64)
if(APPLE)
	set(FMI2_PLATFORM darwin64)
elseif(WIN32)
	set(FMI2_PLATFORM win64)
endif()
set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/fmi2_linear")
file(MAKE_DIRECTORY "${STAGE}")
configure_file(modelDescription2.xml "${STAGE}/modelDescription.xml" COPYONLY)
add_library(fmi_fixture_fmi2 MODULE TestModel2.cpp)
target_include_directories(fmi_fixture_fmi2 PRIVATE ../third_party/fmi2)
set_target_properties(fmi_fixture_fmi2 PROPERTIES
	PREFIX ""
	OUTPUT_NAME TestModel2
	SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI2_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI2_PLATFORM}"
)
set(FMU "${CMAKE_CURRENT_BINARY_DIR}/fmi2_linear.fmu")
add_custom_command(OUTPUT "${FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_fmi2 "${STAGE}/modelDescription.xml"
)
add_custom_target(fmi_package_fmi2_linear DEPENDS "${FMU}")
add_dependencies(FMI_ImportTest fmi_package_fmi2_linear)
add_test(NAME fmi_fmi2_linear COMMAND FMI_ImportTest "${FMU}" fmi2_linear)

set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/fmi2_nonlinear_threephase")
file(MAKE_DIRECTORY "${STAGE}")
configure_file(modelDescription2ThreePh.xml "${STAGE}/modelDescription.xml" COPYONLY)
add_library(fmi_fixture_fmi2_nonlinear_threephase MODULE TestModel2ThreePh.cpp)
target_include_directories(fmi_fixture_fmi2_nonlinear_threephase PRIVATE ../third_party/fmi2)
set_target_properties(fmi_fixture_fmi2_nonlinear_threephase PROPERTIES
	PREFIX ""
	OUTPUT_NAME TestModel2ThreePh
	SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI2_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI2_PLATFORM}"
)
set(FMU "${CMAKE_CURRENT_BINARY_DIR}/fmi2_nonlinear_threephase.fmu")
add_custom_command(OUTPUT "${FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_fmi2_nonlinear_threephase "${STAGE}/modelDescription.xml"
)
add_custom_target(fmi_package_fmi2_nonlinear_threephase DEPENDS "${FMU}")
add_executable(FMI_NonlinearTest2 NonlinearTest.cpp)
target_link_libraries(FMI_NonlinearTest2 PRIVATE dpsim-fmi)
target_include_directories(FMI_NonlinearTest2 PRIVATE ../examples/common)
add_dependencies(FMI_NonlinearTest2 fmi_package_fmi2_nonlinear_threephase)
add_test(NAME fmi_fmi2_nonlinear_threephase
	COMMAND FMI_NonlinearTest2 "${FMU}")

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
target_include_directories(FMI_ImportTest3Ph PRIVATE ../examples/common)
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
target_include_directories(FMI_NonlinearTest PRIVATE ../examples/common)
add_dependencies(FMI_NonlinearTest fmi_package_nonlinear_threephase)

add_test(NAME fmi_nonlinear_threephase COMMAND FMI_NonlinearTest "${FMU}")
add_test(NAME fmi_nonlinear_reuse_jacobian COMMAND FMI_NonlinearTest "${FMU}" reuse)
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
target_include_directories(FMI_EventCurrentTest PRIVATE ../examples/common)
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
target_include_directories(FMI_StateEventTest PRIVATE ../examples/common)
add_dependencies(FMI_StateEventTest fmi_package_state_event)

add_test(NAME fmi_state_event_reset COMMAND FMI_StateEventTest "${FMU}")

set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/sampled_threephase")
file(MAKE_DIRECTORY "${STAGE}")

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/modelDescription3Ph.xml" SAMPLED_XML)

string(REPLACE "name=\"u_a\" valueReference=\"0\" causality=\"input\" variability=\"continuous\" start=\"0\""
	"name=\"u_a\" valueReference=\"0\" causality=\"input\" variability=\"continuous\" start=\"7\"" SAMPLED_XML "${SAMPLED_XML}")
string(REPLACE "name=\"u_b\" valueReference=\"1\" causality=\"input\" variability=\"continuous\" start=\"0\""
	"name=\"u_b\" valueReference=\"1\" causality=\"input\" variability=\"continuous\" start=\"-11\"" SAMPLED_XML "${SAMPLED_XML}")
string(REPLACE "name=\"u_c\" valueReference=\"2\" causality=\"input\" variability=\"continuous\" start=\"0\""
	"name=\"u_c\" valueReference=\"2\" causality=\"input\" variability=\"continuous\" start=\"13\"" SAMPLED_XML "${SAMPLED_XML}")

string(
	REPLACE
		"<ModelExchange "
		"<ModelExchange canGetAndSetFMUState=\"true\" "
		SAMPLED_XML
		"${SAMPLED_XML}"
)

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
target_compile_definitions(fmi_fixture_sampled_threephase PRIVATE FIXTURE_NONLINEAR=1 FIXTURE_SAMPLED=1 FIXTURE_SNAPSHOTS=1)

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
target_include_directories(FMI_SampledTest PRIVATE ../examples/common)
add_dependencies(FMI_SampledTest fmi_package_sampled_threephase)

add_test(NAME fmi_sampled_threephase COMMAND FMI_SampledTest "${FMU}")
add_test(NAME fmi_sampled_reuse_jacobian COMMAND FMI_SampledTest "${FMU}" reuse)
add_test(NAME fmi_multirate_sampled_tick COMMAND FMI_SampledTest "${FMU}" multirate)

# A second instance has a 0.5 ms event clock, proving that each FMU keeps its
# own deterministic schedule while both are synchronized to one EMT macro step.
set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/sampled_threephase_fast")
file(MAKE_DIRECTORY "${STAGE}")
file(WRITE "${STAGE}/modelDescription.xml" "${SAMPLED_XML}")
add_library(fmi_fixture_sampled_threephase_fast MODULE TestModel3Ph.cpp)
target_include_directories(fmi_fixture_sampled_threephase_fast PRIVATE ../third_party/fmi3)
target_compile_definitions(fmi_fixture_sampled_threephase_fast PRIVATE
	FIXTURE_NONLINEAR=1 FIXTURE_SAMPLED=1 FIXTURE_SNAPSHOTS=1
	FIXTURE_SAMPLED_PERIOD=0.0005)
set_target_properties(fmi_fixture_sampled_threephase_fast PROPERTIES
	PREFIX "" OUTPUT_NAME TestModel3Ph SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
)
set(FAST_FMU "${CMAKE_CURRENT_BINARY_DIR}/sampled_threephase_fast.fmu")
add_custom_command(OUTPUT "${FAST_FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${FAST_FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_sampled_threephase_fast "${STAGE}/modelDescription.xml"
)
add_custom_target(fmi_package_sampled_threephase_fast DEPENDS "${FAST_FMU}")
add_executable(FMI_MultiRateTest MultiRateTest.cpp)
target_link_libraries(FMI_MultiRateTest PRIVATE dpsim-fmi)
target_include_directories(FMI_MultiRateTest PRIVATE ../examples/common)
add_dependencies(FMI_MultiRateTest
	fmi_package_sampled_threephase fmi_package_sampled_threephase_fast)
add_test(NAME fmi_multirate_distinct_fmu_clocks
	COMMAND FMI_MultiRateTest "${CMAKE_CURRENT_BINARY_DIR}/sampled_threephase.fmu" "${FAST_FMU}")
add_test(NAME fmi_reject_sampled_lti COMMAND FMI_ImportTest3Ph "${FMU}" reject_sampled)

# The same sampled model without FMI rollback must fail before a crossed event
# is probed or committed.
set(STAGE "${CMAKE_CURRENT_BINARY_DIR}/stage/sampled_no_rollback")
file(MAKE_DIRECTORY "${STAGE}")
string(REPLACE " canGetAndSetFMUState=\"true\"" ""
	SAMPLED_NO_ROLLBACK_XML "${SAMPLED_XML}")
file(WRITE "${STAGE}/modelDescription.xml" "${SAMPLED_NO_ROLLBACK_XML}")

add_library(fmi_fixture_sampled_no_rollback MODULE TestModel3Ph.cpp)
target_include_directories(fmi_fixture_sampled_no_rollback PRIVATE ../third_party/fmi3)
target_compile_definitions(fmi_fixture_sampled_no_rollback PRIVATE FIXTURE_NONLINEAR=1 FIXTURE_SAMPLED=1)
set_target_properties(fmi_fixture_sampled_no_rollback PROPERTIES
	PREFIX "" OUTPUT_NAME TestModel3Ph SUFFIX "${FMI_SUFFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
	RUNTIME_OUTPUT_DIRECTORY "${STAGE}/binaries/${FMI_PLATFORM}"
)
set(NO_ROLLBACK_FMU "${CMAKE_CURRENT_BINARY_DIR}/sampled_no_rollback.fmu")
add_custom_command(OUTPUT "${NO_ROLLBACK_FMU}"
	COMMAND ${CMAKE_COMMAND} -E tar cf "${NO_ROLLBACK_FMU}" --format=zip -- modelDescription.xml binaries
	WORKING_DIRECTORY "${STAGE}"
	DEPENDS fmi_fixture_sampled_no_rollback "${STAGE}/modelDescription.xml"
)
add_custom_target(fmi_package_sampled_no_rollback DEPENDS "${NO_ROLLBACK_FMU}")
add_dependencies(FMI_SampledTest fmi_package_sampled_no_rollback)
add_test(NAME fmi_multirate_requires_rollback
	COMMAND FMI_SampledTest "${NO_ROLLBACK_FMU}" no_rollback)
add_test(NAME fmi_multirate_single_pass
	COMMAND FMI_SampledTest "${NO_ROLLBACK_FMU}" single_pass)

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
