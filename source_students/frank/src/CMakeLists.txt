ADD_EXECUTABLE( noise_generator noise_generator.cpp )
TARGET_LINK_LIBRARIES( noise_generator myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

#find_library( PNG_LIBRARY NAMES png PATHS ENV LD_LIBRARY_PATH DEFAULT_PATH )

#ADD_EXECUTABLE( mm2png mm2png.cpp pngwriter.cc)
#TARGET_LINK_LIBRARIES( mm2png myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} ${PNG_LIBRARY} )

ADD_EXECUTABLE( noise_solver noise_solver.cpp )
TARGET_LINK_LIBRARIES( noise_solver myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

ADD_EXECUTABLE( psd psd.cpp )
TARGET_LINK_LIBRARIES( psd myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

ADD_EXECUTABLE( reduce_data reduce_data.cpp )
TARGET_LINK_LIBRARIES( reduce_data myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

#ADD_EXECUTABLE( lis_sg_test_3 lis_sg_test_3.cpp )
#TARGET_LINK_LIBRARIES( lis_sg_test_3 myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

ADD_EXECUTABLE( noise_int_demo noise_int_demo.cpp )
TARGET_LINK_LIBRARIES( noise_int_demo myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

ADD_EXECUTABLE( noise_int_demo_B noise_int_demo_B.cpp )
TARGET_LINK_LIBRARIES( noise_int_demo_B myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

ADD_EXECUTABLE( prob_dist prob_dist.cpp )
TARGET_LINK_LIBRARIES( prob_dist myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )

ADD_EXECUTABLE( noise_gen noise_gen.cpp )
TARGET_LINK_LIBRARIES( noise_gen myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${MUPARSER_LIBRARY} ${LIS_LIBRARY} )
