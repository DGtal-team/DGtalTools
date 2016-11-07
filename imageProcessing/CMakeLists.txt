
# These tools require Eigen (because DEC requires Eigen)
if(WITH_EIGEN)
  SET(DGTAL_TOOLS_SRC
    at-u0-v1
    at-u2-v0
    )
  FOREACH(FILE ${DGTAL_TOOLS_SRC})
    add_executable(${FILE} ${FILE})
    target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} )
    install(TARGETS ${FILE}
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
  ENDFOREACH(FILE)
endif(WITH_EIGEN)