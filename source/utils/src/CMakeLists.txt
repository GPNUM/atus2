
# ADD_EXECUTABLE( Parameterberechnung Parameterberechnung.cpp)
# TARGET_LINK_LIBRARIES( Parameterberechnung myutils)

ADD_EXECUTABLE( ana_tool ana_tool.cpp )
TARGET_LINK_LIBRARIES( ana_tool myutils ${MUPARSER_LIBRARY} )

ADD_EXECUTABLE( ana_expval ana_expval.cpp )
TARGET_LINK_LIBRARIES( ana_expval myutils )

ADD_EXECUTABLE( ana_expval_ex ana_expval_ex.cpp )
TARGET_LINK_LIBRARIES( ana_expval_ex myutils )

ADD_EXECUTABLE( ana_overlap ana_overlap.cpp )
TARGET_LINK_LIBRARIES( ana_overlap myutils )

ADD_EXECUTABLE( gen_mirror gen_mirror.cpp )
TARGET_LINK_LIBRARIES( gen_mirror myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m )

ADD_EXECUTABLE( ana_impuls ana_impulsverteilung.cpp )
TARGET_LINK_LIBRARIES( ana_impuls myutils )

ADD_EXECUTABLE( gen_psi_0 gen_psi_0.cpp )
TARGET_LINK_LIBRARIES( gen_psi_0 myutils ${MUPARSER_LIBRARY} )

ADD_EXECUTABLE( inflate_domain inflate_domain.cpp )
TARGET_LINK_LIBRARIES( inflate_domain myutils )

ADD_EXECUTABLE( compute_overlap compute_overlap.cpp)
TARGET_LINK_LIBRARIES( compute_overlap myutils)

#ADD_EXECUTABLE( compute_expansion_coefficients compute_expansion_coefficients.cpp )
#TARGET_LINK_LIBRARIES( compute_expansion_coefficients myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} )

ADD_EXECUTABLE( gen_head gen_head.cpp )
TARGET_LINK_LIBRARIES( gen_head myutils )

ADD_EXECUTABLE(gen_pvd gen_pvd.cpp)
