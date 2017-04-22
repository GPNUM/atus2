
### atom interferometer 1 species
ADD_EXECUTABLE( bragg bragg.cpp  )
TARGET_LINK_LIBRARIES( bragg myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( double_bragg double_bragg.cpp  )
TARGET_LINK_LIBRARIES( double_bragg myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( raman raman.cpp  )
TARGET_LINK_LIBRARIES( raman myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( bragg_ds bragg_ds.cpp  )
TARGET_LINK_LIBRARIES( bragg_ds myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})
