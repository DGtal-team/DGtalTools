
SET(DGTAL_TOOLS_SRC
  pgm2freeman
  freeman2sdp
  vol2raw
  raw2vol
  ofs2off
  convertVol
  vol2sdp
  slice2vol
  vol2slice
  vol2off
  vol2obj
  longvol2vol)

if( WITH_HDF5 )
SET(DGTAL_TOOLS_SRC ${DGTAL_TOOLS_SRC} 
    raw2HDF5
    HDF52vol)
endif( WITH_HDF5 )

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} )
  install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)


if (  WITH_ITK )
  SET(DGTAL_TOOLS_ITK
    dicom2vol)
  FOREACH(FILE ${DGTAL_TOOLS_ITK})
     add_executable(${FILE} ${FILE})
     target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} )
     install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)
ENDIF( WITH_ITK )
