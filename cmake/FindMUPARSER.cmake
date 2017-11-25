
# Find the MUPARSER library
# Usage:
#   find_package( MUPARSER [REQUIRED] [QUIET] )
#

if( MUPARSER_ROOT )

    find_library(  MUPARSER_LIBRARY
                   NAMES "muparser" 
                   PATHS ${MUPARSER_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_path(  MUPARSER_INCLUDE_DIR
                NAMES "muPaser.h"
                PATHS ${MUPARSER_ROOT}
                PATH_SUFFIXES "include"
                NO_DEFAULT_PATH
    )

else()

    find_library(   MUPARSER_LIBRARY
                    NAMES "muparser"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    get_filename_component( TMP ${MUPARSER_LIBRARY} PATH )
    get_filename_component( TMP ${TMP} PATH )
    set( MUPARSER_INCLUDE_DIR ${TMP}/include CACHE STRING INTERNAL )

endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MUPARSER
      REQUIRED_VARS MUPARSER_INCLUDE_DIR MUPARSER_LIBRARY
      HANDLE_COMPONENTS
      )

mark_as_advanced(
        MUPARSER_LIBRARY
        MUPARSER_INCLUDE_DIR
      )
