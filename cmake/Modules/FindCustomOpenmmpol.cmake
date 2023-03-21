# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2023  DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomOpenmmpol
----------------

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
  the settings in ``OPENMMPOL_LIBRARY`` will be used without any further
  checks.

  ``OPENMMPOL_LIBRARY``

  Customized openmmpol library/libraries to use (instead of autodetected one). If
  not set or empty, the default library (openmmpol) will be tried. The listed
  libraries will be checked for existence (unless disabled in
  ``OPENMMPOL_DETECTION``) and the variable is overwritten to contain the libraries
  with their with full path.

  ``OPENMMPOL_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.


#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET Openmmpol::Openmmpol)

  set(CUSTOM_OPENMMPOL_FOUND True)
  set(CustomOpenmmpol_FOUND True)
  set(OPENMMPOL_FOUND True)
  set(Openmmpol_FOUND True)

else()

  option(OPENMMPOL_DETECTION "Whether openmmpol library should be detected" TRUE)

  if(OPENMMPOL_DETECTION)

    find_package(PkgConfig)
    pkg_check_modules(_openmmpol QUIET openmmpol)

    # Overwrite PkgConfig values by user defined input if present.
    if(NOT "${OPENMMPOL_LIBRARY}" STREQUAL "")
      set(_openmmpol_LIBRARIES ${OPENMMPOL_LIBRARY})
      set(_openmmpol_LIBRARY_DIRS ${OPENMMPOL_LIBRARY_DIR})
    endif()

    find_custom_libraries("${_openmmpol_LIBRARIES}" "${_openmmpol_LIBRARY_DIRS}"
      "${CustomOpenmmpol_FIND_QUIETLY}" _libs)
    set(OPENMMPOL_LIBRARY "${_libs}" CACHE STRING "List of openmmpol libraries to link" FORCE)
    unset(_libs)
    unset(_openmmpol_LIBRARIES)
    unset(_openmmpol_LIBRARY_DIRS)

    set(OPENMMPOL_DETECTION False CACHE BOOL "Whether openmmpol libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomOpenmmpol REQUIRED_VARS OPENMMPOL_LIBRARY)

  # find_library(CUSTOM_OPENMMPOL NAMES libopenmmpol.so HINTS ${OPENMMPOL_LIBRARY_DIR}) 
  # add_library(Openmmpol::Openmmpol INTERFACE IMPORTED)
  # set_target_properties(Openmmpol::Openmmpol PROPERTIES IMPORTED_LOCATION ${CUSTOM_OPENMMPOL})

  set(OPENMMPOL_FOUND ${CUSTOMOPENMMPOL_FOUND})
  set(Openmmpol_FOUND ${CUSTOMOPENMMPOL_FOUND})

  if(OPENMMPOL_FOUND AND NOT TARGET Openmmpol::Openmmpol)
    add_library(Openmmpol::Openmmpol INTERFACE IMPORTED)
    target_link_libraries(Openmmpol::Openmmpol INTERFACE "${OPENMMPOL_LIBRARY}")
  endif()

endif()
