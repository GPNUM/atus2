
ADD_EXECUTABLE( bragg_mpi bragg_mpi.cpp)
TARGET_LINK_LIBRARIES( bragg_mpi myutils_mpi_2 myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( double_bragg_mpi double_bragg_mpi.cpp)
TARGET_LINK_LIBRARIES( double_bragg_mpi myutils_mpi_2 myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( bragg_ds_mpi bragg_ds_mpi.cpp)
TARGET_LINK_LIBRARIES( bragg_ds_mpi myutils_mpi_2 myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( double_bragg_ds_mpi double_bragg_ds_mpi.cpp)
TARGET_LINK_LIBRARIES( double_bragg_ds_mpi myutils_mpi_2 myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})
