
ADD_EXECUTABLE( rt_solver_mpi rt_solver_mpi.cpp )
TARGET_LINK_LIBRARIES( rt_solver_mpi myutils_mpi_2 myutils ${MUPARSER_LIBRARY} )

ADD_EXECUTABLE( sobmin_mpi sob_solver_mpi.cpp  )
TARGET_LINK_LIBRARIES( sobmin_mpi myutils_mpi_2 myutils ${MUPARSER_LIBRARY} )

