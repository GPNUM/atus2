
### real time propagation 
ADD_EXECUTABLE( rt_solver rt_solver.cpp  )
TARGET_LINK_LIBRARIES( rt_solver myutils ${MUPARSER_LIBRARY} )

### STATIONARY SOLVER
ADD_EXECUTABLE( sobmin sob_solver.cpp )
TARGET_LINK_LIBRARIES( sobmin myutils ${MUPARSER_LIBRARY} )

ADD_EXECUTABLE( sobmin_2 sob_solver_2.cpp )
TARGET_LINK_LIBRARIES( sobmin_2 myutils ${MUPARSER_LIBRARY} )
