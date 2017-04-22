
ADD_LIBRARY( myutils_mpi_2 cft_2d_MPI.cpp cft_3d_MPI.cpp timer.cpp )
TARGET_LINK_LIBRARIES( myutils_mpi_2 m ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARY_3} ${FFTW_LIBRARY_2} ${FFTW_LIBRARY_1}  )

ADD_EXECUTABLE( fftw_3d_mpi_bench fftw3_3d_MPI_bench.cpp )
TARGET_LINK_LIBRARIES( fftw_3d_mpi_bench myutils_mpi_2 )

ADD_EXECUTABLE( mpi_test MPI_test.cpp )
TARGET_LINK_LIBRARIES( mpi_test myutils_mpi_2 )
  
ADD_EXECUTABLE( fftw_2d_mpi_test fftw3_2d_MPI_test.cpp )
TARGET_LINK_LIBRARIES( fftw_2d_mpi_test myutils_mpi_2 )

ADD_EXECUTABLE( fftw_3d_mpi_test fftw3_3d_MPI_test.cpp )
TARGET_LINK_LIBRARIES( fftw_3d_mpi_test myutils_mpi_2 )

