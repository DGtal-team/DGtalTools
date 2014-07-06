SET(DGTAL_TOOLS_SRC
	displayContours
        patternTriangulation
)

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
  install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)

if (  WITH_VISU3D_QGLVIEWER )

  SET(DGTAL_EXAMPLES_QGL_SRC

    3dVolViewer
    3dVolBoundaryViewer
    3dCurvatureViewer
    3dCurvatureViewerNoise
    3dCurveViewer
    3dSDPViewer
    meshViewer)
  FOREACH(FILE ${DGTAL_EXAMPLES_QGL_SRC}) 
    add_executable(${FILE} ${FILE})
    target_link_libraries ( ${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})   
   
    install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
  ENDFOREACH(FILE)
  add_executable(3dImageViewer 3dImageViewer specificClasses/Viewer3DImage  )
  target_link_libraries ( 3dImageViewer  ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})   
  install(TARGETS 3dImageViewer 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
QT4_WRAP_UI(ui_2d3dImageViewer_FORMS_HEADERS 2d3dImageViewer.ui)
QT4_WRAP_CPP(2d3dImageViewer_HEADERS_MOC 2d3dImageViewer.h OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)


ADD_EXECUTABLE(2d3dImageViewer 2d3dImageViewer.cpp  ${2d3dImageViewer_HEADERS_MOC}  
    ${ui_2d3dImageViewer_FORMS_HEADERS} )
 target_link_libraries ( 2d3dImageViewer  ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})   

endif ( WITH_VISU3D_QGLVIEWER   )



