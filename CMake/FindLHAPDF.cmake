# Try to find LHAPDF
# Defines:
#   LHAPDF_FOUND
#   LHAPDF_INCLUDE_DIR
#   LHAPDF_INCLUDE_DIRS (not cached)
#   LHAPDF_LIBRARY
#   LHAPDF_LIBRARIES (not cached)
#   LHAPDF_LIBRARY_DIR (not cached)

find_library(LHAPDF_LIBRARY NAMES LHAPDF
    HINTS $ENV{LHAPDF_ROOT_DIR}/lib ${LHAPDF_ROOT_DIR}/lib)

IF(${LHAPDF_LIBRARY} MATCHES "LHAPDF_LIBRARY-NOTFOUND")
    FIND_PROGRAM(LHAPDF_CONFIG_EXECUTABLE NAMES lhapdf-config
        HINTS $ENV{LHAPDF_ROOT_DIR}/bin ${LHAPDF_ROOT_DIR}/bin $ENV{PATH})
    IF(${LHAPDF_CONFIG_EXECUTABLE} MATCHES "LHAPDF_CONFIG_EXECUTABLE-NOTFOUND")
        MESSAGE(STATUS "Looking for LHAPDF... lhapdf-config executable not found")
    ELSE(${LHAPDF_CONFIG_EXECUTABLE} MATCHES "LHAPDF_CONFIG_EXECUTABLE-NOTFOUND")
        MESSAGE(STATUS "Looking for LHAPDF... using lhapdf-config executable")
        EXEC_PROGRAM(${LHAPDF_CONFIG_EXECUTABLE} ARGS "--prefix" OUTPUT_VARIABLE LHAPDF_PREFIX)
        find_library(LHAPDF_LIBRARY NAMES "LHAPDF" PATHS ${LHAPDF_PREFIX}/lib)
    ENDIF(${LHAPDF_CONFIG_EXECUTABLE} MATCHES "LHAPDF_CONFIG_EXECUTABLE-NOTFOUND")
ENDIF(${LHAPDF_LIBRARY} MATCHES "LHAPDF_LIBRARY-NOTFOUND")

find_path(LHAPDF_INCLUDE_DIR LHAPDF/LHAPDF.h
    HINTS $ENV{LHAPDF_ROOT_DIR}/include ${LHAPDF_ROOT_DIR}/include ${LHAPDF_PREFIX}/include)

mark_as_advanced(LHAPDF_LIBRARY LHAPDF_INCLUDE_DIR)

# handle QUIETLY and REQUIRED arguments and set LHAPDF_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LHAPDF DEFAULT_MSG LHAPDF_INCLUDE_DIR LHAPDF_LIBRARY)

set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARY})
get_filename_component(LHAPDF_LIBRARY_DIRS ${LHAPDF_LIBRARY} PATH)

set(LHAPDF_INCLUDE_DIRS ${LHAPDF_INCLUDE_DIR})

mark_as_advanced(LHAPDF_FOUND)