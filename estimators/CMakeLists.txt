SET(DGTAL_TOOLS_SRC
  2dLocalEstimators
  lengthEstimators
  vol2normalField
  curvatureBC
  curvatureMCMS
  tangentBC
  statisticsEstimators
  curvatureScaleSpaceBCC
)

if (  WITH_CGAL )
  SET(CGAL_TESTS_SRC
    3dLocalEstimators )
  FOREACH(FILE ${CGAL_TESTS_SRC})
    add_executable(${FILE} ${FILE})
    target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
  ENDFOREACH(FILE)
endif ( WITH_CGAL )

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
  install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)
