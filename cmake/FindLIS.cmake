
# Find the LIS library http://www.ssisc.org/lis/index.en.html
# Usage:
#   find_package( LIS [REQUIRED] [QUIET] )
#

if( LIS_ROOT )

    find_library(  LIS_LIBRARY
                   NAMES "lis" 
                   PATHS ${LIS_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_path(  LIS_INCLUDE_DIR
                NAMES "lis.h"
                PATHS ${LIS_ROOT}
                PATH_SUFFIXES "include"
                NO_DEFAULT_PATH
    )

else()

    find_library(   LIS_LIBRARY
                    NAMES "lis"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    get_filename_component( TMP ${LIS_LIBRARY} PATH )
    get_filename_component( TMP ${TMP} PATH )
    set( LIS_INCLUDE_DIR ${TMP}/include CACHE STRING INTERNAL )

endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LIS
      REQUIRED_VARS LIS_INCLUDE_DIR LIS_LIBRARY
      HANDLE_COMPONENTS
      )

mark_as_advanced(
      LIS_LIBRARY
      LIS_INCLUDE_DIR
      )
