SET(DGTAL_TOOLS_SRC
  2dLocalEstimators
  lengthEstimators
  vol2normalField
  curvatureBC
  curvatureMCMS
  tangentBC
  statisticsEstimators
  generic3dNormalEstimators
  curvatureScaleSpaceBCC
  eulerCharacteristic
)

if (WITH_EIGEN)
    add_executable(volSurfaceRegularization volSurfaceRegularization.cpp volSurfaceRegularization-details/surface_approx.cpp)
    target_link_libraries(volSurfaceRegularization ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
endif()

if ( WITH_VISU3D_QGLVIEWER )
  SET(VISU3D_TESTS_SRC
      3dCurveTangentEstimator
  )
  FOREACH(FILE ${VISU3D_TESTS_SRC})
    add_executable(${FILE} ${FILE})
    target_link_libraries (${FILE}  ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
  ENDFOREACH(FILE)
endif ( WITH_VISU3D_QGLVIEWER )

if (  WITH_CGAL )
  SET(CGAL_TESTS_SRC
#    3dLocalEstimators
     )
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
