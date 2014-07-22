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
    meshViewer
    3dHeightMapViewer
    3dCompSurfelData)

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
QT4_WRAP_UI(ui_sliceViewer_FORMS_HEADERS sliceViewer.ui)
QT4_WRAP_CPP(sliceViewer_HEADERS_MOC sliceViewer.h OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)


ADD_EXECUTABLE(sliceViewer sliceViewer.cpp  ${sliceViewer_HEADERS_MOC}  
    ${ui_sliceViewer_FORMS_HEADERS} )
 target_link_libraries ( sliceViewer  ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})   

endif ( WITH_VISU3D_QGLVIEWER   )



