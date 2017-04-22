

ADD_EXECUTABLE( batch_case_generation_4_tammy batch_case_generation_4_tammy.cpp)
TARGET_LINK_LIBRARIES( batch_case_generation_4_tammy myutils)

ADD_EXECUTABLE( compute_rabi_times compute_rabi_times.cpp)
TARGET_LINK_LIBRARIES(compute_rabi_times myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( tammy_auswertung compute_rabi_times_4_tammy.cpp)
TARGET_LINK_LIBRARIES(tammy_auswertung myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( tammy_diff_auswertung compute_rabi_times_diff_4_tammy.cpp)
TARGET_LINK_LIBRARIES(tammy_diff_auswertung myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_EXECUTABLE( compute_expval_ex_diff_4_tammy compute_expval_ex_diff_4_tammy.cpp)
TARGET_LINK_LIBRARIES(compute_expval_ex_diff_4_tammy myutils ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})
