# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2023  DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomOpenmmpol
---------------

Finds the openmmpol library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``Openmmpol::Openmmpol``
  The openmmpol library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``OPENMMPOL_FOUND``
  True if the system has the openmmpol library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``OPENMMPOL_DETECTION``
  Whether openmmpol libraries should be detected (default: True). If set to False,
  the settings in ``OPENMMPOL_LIBRARY`` and ``OPENMMPOL_INCLUDE_DIR`` will be used without any further
  checks.

``OPENMMPOL_LIBRARY``

  Customized openmmpol library/libraries to use (instead of autodetected one). If
  not set or empty, the default library (openMMPol) will be tried. The listed
  libraries will be checked for existence (unless disabled in
  ``OPENMMPOL_DETECTION``) and the variable is overwritten to contain the libraries
  with their with full path.

``OPENMMPOL_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.

``OPENMMPOL_INCLUDE_DIRECTORY``
  Customized openmmpol include directory.

#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET Openmmpol::Openmmpol)

  set(CUSTOMOPENMMPOL_FOUND True)
  set(CustomOpenmmpol_FOUND True)
  set(OPENMMPOL_FOUND True)
  set(Openmmpol_FOUND True)

else()

  option(OPENMMPOL_DETECTION "Whether openmmpol library should be detected" TRUE)

  if(OPENMMPOL_DETECTION)

    if("${OPENMMPOL_LIBRARY}" STREQUAL "")

      if (NOT OPENMMPOL_FOUND)
        find_package(openmmpol 0.2 QUIET)
      endif()

      set(OPENMMPOL_LIBRARY "${OPENMMPOL_LIBRARIES}" CACHE STRING "openmmpol library to link" FORCE)
      set(OPENMMPOL_LINKER_FLAG "${OPENMMPOL_LINKER_FLAGS}" CACHE STRING 
          "Linker flags to use when linking openmmpol" FORCE)
      set(OPENMMPOL_INCLUDE_DIRECTORY "${OPENMMPOL_INCLUDE_DIRS}" CACHE STRING "openmmpol include directory" FORCE)

    elseif(NOT "${OPENMMPOL_LIBRARY}" STREQUAL "NONE")

      find_custom_libraries("${OPENMMPOL_LIBRARY}" "${OPENMMPOL_LIBRARY_DIR}"
        "${CustomOpenmmpol_FIND_QUIETLY}" _libs)
      set(OPENMMPOL_LIBRARY "${_libs}" CACHE STRING "List of openmmpol libraries to link" FORCE)
      unset(_libs)

    endif()

    set(OPENMMPOL_DETECTION False CACHE BOOL "Whether openmmpol libraries should be detected" FORCE)

    #find_package(PkgConfig QUIET)
    #pkg_check_modules(_openmmpol REQUIRED openmmpol)

    # cmake_print_variables(openmmpol_LIBRARY)
    # # Overwrite PkgConfig values by user defined input if present.
    # if(NOT "${OPENMMPOL_LIBRARY}" STREQUAL "")
    #   set(_openmmpol_LIBRARIES ${OPENMMPOL_LIBRARY})
    #   set(_openmmpol_LIBRARY_DIRS ${OPENMMPOL_LIBRARY_DIR})
    # endif()
    # if(NOT "${OPENMMPOL_INCLUDE_DIR}" STREQUAL "")
    #   set(_openmmpol_INCLUDE_DIRS ${OPENMMPOL_INCLUDE_DIR})
    # endif()

    # find_custom_libraries("${_openmmpol_LIBRARIES}" "${_openmmpol_LIBRARY_DIRS}"
    #   "${CustomOpenmmpol_FIND_QUIETLY}" _libs)

    # set(OPENMMPOL_LIBRARY "${_libs}" CACHE STRING "List of openmmpol libraries to link" FORCE)
    # unset(_libs)
    # unset(_openmmpol_LIBRARIES)
    # unset(_openmmpol_LIBRARY_DIRS)

    # # Check include file
    # find_path(OPENMMPOL_INCLUDE_DIRECTORY NAMES ommp_interface.mod PATHS ${_openmmpol_INCLUDE_DIRS}
    #   PATH_SUFFIXES include openMMPol)
    # unset(_openmmpol_INCLUDE_DIRS)

    # set(OPENMMPOL_DETECTION FALSE CACHE BOOL "Whether openmmpol library should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomOpenmmpol REQUIRED_VARS OPENMMPOL_LIBRARY OPENMMPOL_INCLUDE_DIRECTORY)

  set(OPENMMPOL_FOUND ${CUSTOMOPENMMPOL_FOUND})
  set(Openmmpol_FOUND ${CUSTOMOPENMMPOL_FOUND})

  if(OPENMMPOL_FOUND AND NOT TARGET Openmmpol::Openmmpol)
    add_library(Openmmpol::Openmmpol INTERFACE IMPORTED)
    target_link_libraries(Openmmpol::Openmmpol INTERFACE "${OPENMMPOL_LIBRARY}")
    target_include_directories(Openmmpol::Openmmpol INTERFACE "${OPENMMPOL_INCLUDE_DIRECTORY}")

    if(NOT "${OPENMMPOL_LINKER_FLAG}" STREQUAL "")
      target_link_options(Openmmpol::Openmmpol INTERFACE "${OPENMMPOL_LINKER_FLAG}")
    endif()

  endif()

  mark_as_advanced(OPENMMPOL_DETECTION OPENMMPOL_LIBRARY OPENMMPOL_LIBRARY_DIR OPENMMPOL_LINKER_FLAG)
endif()