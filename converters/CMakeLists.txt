
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
  longvol2vol)

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} )
  install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)
