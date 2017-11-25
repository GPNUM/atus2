
# Find the HDF5 library
# Usage:
#   find_package( HDF5 [REQUIRED] [QUIET] )
#

if( HDF5_ROOT )

    find_library(  HDF5_LIBRARY_1
                   NAMES "hdf5_hl_cpp" 
                   PATHS ${HDF5_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_library(  HDF5_LIBRARY_2
                   NAMES "hdf5_cpp"
                   PATHS ${HDF5_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_library(  HDF5_LIBRARY_3
                   NAMES "hdf5_hl"
                   PATHS ${HDF5_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )    
    
    find_library(  HDF5_LIBRARY_4
                   NAMES "hdf5"
                   PATHS ${HDF5_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )    

    find_path(  HDF5_INCLUDE_DIR
                NAMES "hdf5.h"
                PATHS ${HDF5_ROOT}
                PATH_SUFFIXES "include"
                NO_DEFAULT_PATH
    )

else()

    find_library(   HDF5_LIBRARY_1 
                    NAMES "hdf5_hl_cpp"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    find_library(   HDF5_LIBRARY_2 
                    NAMES "hdf5_cpp"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    find_library(   HDF5_LIBRARY_3
                    NAMES "hdf5_hl"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )   
    
    find_library(   HDF5_LIBRARY_4
                    NAMES "hdf5"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )       

    get_filename_component( TMP ${HDF5_LIBRARY_1} PATH )
    get_filename_component( TMP ${TMP} PATH )
    set( HDF5_INCLUDE_DIR ${TMP}/include CACHE STRING INTERNAL )

endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HDF5
      REQUIRED_VARS HDF5_INCLUDE_DIR HDF5_LIBRARY_1 HDF5_LIBRARY_2 HDF5_LIBRARY_3 HDF5_LIBRARY_4
      HANDLE_COMPONENTS
      )

mark_as_advanced(
      HDF5_LIBRARY_1
      HDF5_LIBRARY_2
      HDF5_LIBRARY_3
      HDF5_LIBRARY_4
      HDF5_INCLUDE_DIR
      )
