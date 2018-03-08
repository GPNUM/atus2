
ADD_LIBRARY( myutils cft_1d.cpp cft_2d.cpp cft_3d.cpp rft_1d.cpp rft_2d.cpp rft_3d.cpp misc.cpp noise3_2d.cpp ParameterHandler.cpp zernike.cpp pugixml.cpp )
TARGET_LINK_LIBRARIES( myutils m gomp ${FFTW_LIBRARY_1} ${FFTW_LIBRARY_2} )

ADD_EXECUTABLE( slice_3d slice_3d.cpp )
TARGET_LINK_LIBRARIES( slice_3d myutils )

ADD_EXECUTABLE( test_noise3 test_noise3.cpp )
TARGET_LINK_LIBRARIES( test_noise3 myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m )

ADD_EXECUTABLE( cft_1d_test cft_1d_test.cpp )
TARGET_LINK_LIBRARIES( cft_1d_test myutils m )

ADD_EXECUTABLE( cft_2d_test cft_2d_test.cpp )
TARGET_LINK_LIBRARIES( cft_2d_test myutils m )

ADD_EXECUTABLE( cft_3d_test cft_3d_test.cpp )
TARGET_LINK_LIBRARIES( cft_3d_test myutils m )

ADD_EXECUTABLE( rft_1d_test rft_1d_test.cpp )
TARGET_LINK_LIBRARIES( rft_1d_test myutils m )

ADD_EXECUTABLE( rft_2d_test rft_2d_test.cpp )
TARGET_LINK_LIBRARIES( rft_2d_test myutils m )

ADD_EXECUTABLE( rft_3d_test rft_3d_test.cpp )
TARGET_LINK_LIBRARIES( rft_3d_test myutils m )
