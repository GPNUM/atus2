
# Find the GSL library
# Usage:
#   find_package( GSL [REQUIRED] [QUIET] )
#

if( GSL_ROOT )

    find_library(  GSL_LIBRARY_1
                   NAMES "gsl" 
                   PATHS ${GSL_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_library(  GSL_LIBRARY_2
                   NAMES "gslcblas"
                   PATHS ${GSL_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_path(  GSL_INCLUDE_DIR
                NAMES "gsl_version.h"
                PATHS ${GSL_ROOT}
                PATH_SUFFIXES "include"
                NO_DEFAULT_PATH
    )

else()

    find_library(   GSL_LIBRARY_1 
                    NAMES "gsl"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    find_library(   GSL_LIBRARY_2 
                    NAMES "gslcblas"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    get_filename_component( TMP ${GSL_LIBRARY_1} PATH )
    get_filename_component( TMP ${TMP} PATH )
    set( GSL_INCLUDE_DIR ${TMP}/include CACHE STRING INTERNAL )

endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GSL
      REQUIRED_VARS GSL_INCLUDE_DIR GSL_LIBRARY_1 GSL_LIBRARY_2 
      HANDLE_COMPONENTS
      )

mark_as_advanced(
      GSL_LIBRARY_1
      GSL_LIBRARY_2
      GSL_INCLUDE_DIR
      )
