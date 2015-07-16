#
#  Try to find CFITSIO library
#
#  The module defines:
#       CFITSIO_FOUND        - system has CFITSIO library
#       CFITSIO_LIBRARIES    - libraries (cfitsio and m (if UNIX))
#       CFITSIO_INCLUDE_DIR  - path to header file
#

include(FindPackageHandleStandardArgs)

if (UNIX) # assume on the UNIX platform the CFITSIO library was installed in standard way
    find_path (CFITSIO_INCLUDE_DIR fitsio.h)
    find_library (CFITSIO_LIBRARY NAMES cfitsio PATH_SUFFIXES lib lib64)
    # CFITSIO needs -lm
    find_library(MATH_LIBRARY m)
    find_package_handle_standard_args(CFITSIO  DEFAULT_MSG
                                      CFITSIO_LIBRARY MATH_LIBRARY CFITSIO_INCLUDE_DIR)
    mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY MATH_LIBRARY)
    set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY} ${MATH_LIBRARY})
elseif (WIN32) # user may set CFITSIO_INSTALL_DIR variable to detect header and library files
    set(CFITSIO_INSTALL_DIR "" CACHE STRING "CFITSIO install dir")
    set(CFITSIO_INSTALL_DIR_INTERNAL "" CACHE STRING "CFITSIO install dir")
#    set(CFITSIO_INSTALL_DIR_INTERNAL "g:/PROGRAMS/LIBS/CFITSIO/" CACHE STRING "CFITSIO install dir")
    if(NOT "${CFITSIO_INSTALL_DIR}" STREQUAL "${CFITSIO_INSTALL_DIR_INTERNAL}") # CFITSIO_INSTALL_DIR is given in command-line
        unset(CFITSIO_LIBRARY CACHE)
        unset(CFITSIO_INCLUDE_DIR CACHE)
    endif()
    find_path (CFITSIO_INCLUDE_DIR NAMES fitsio.h PATHS ${CFITSIO_INSTALL_DIR})
    find_library (CFITSIO_LIBRARY NAMES cfitsio PATHS ${CFITSIO_INSTALL_DIR})
    find_package_handle_standard_args(CFITSIO  DEFAULT_MSG
                                      CFITSIO_LIBRARY CFITSIO_INCLUDE_DIR)
    mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY)
    set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY})
endif ()
