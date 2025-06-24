if (CMAKE_VERSION VERSION_LESS 3.9)
  message (FATAL_ERROR "CMake >= 3.9 required")
endif (CMAKE_VERSION VERSION_LESS 3.9)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was suitesparse-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include (CMakeFindDependencyMacro)





# Add the targets after the dependecies were found since SuiteSparse targets
# depend on the external targets.
include ("${CMAKE_CURRENT_LIST_DIR}/suitesparse-targets.cmake")

set (_SuiteSparse_NAMESPACE "")

foreach (_comp ${SuiteSparse_FIND_COMPONENTS})
  set (_TARGET ${_SuiteSparse_NAMESPACE}${_comp})

  if (TARGET ${_TARGET})
    set (SuiteSparse_${_comp}_FOUND TRUE)
  else (TARGET ${_TARGET})
    set (SuiteSparse_${_comp}_FOUND FALSE)
    set (SuiteSparse_NOT_FOUND_MESSAGE "SuiteSparse could not be found because the component ${_comp} could not be found.")
  endif (TARGET ${_TARGET})
endforeach (_comp)

check_required_components (SuiteSparse)
