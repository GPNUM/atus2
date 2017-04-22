
ADD_EXECUTABLE( gpo3 gpo3.cpp )
TARGET_LINK_LIBRARIES( gpo3 myutils )

if( VTK_FOUND )
  ADD_EXECUTABLE( gen_vti gen_vti.cpp )
  TARGET_LINK_LIBRARIES( gen_vti myutils ${VTK_LIBRARIES} )
endif()