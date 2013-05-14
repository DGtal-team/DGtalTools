SET(DGTAL_TOOLS_SRC
	displayContours
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
    3dCurvatureViewer
    3dCurveViewer
    meshViewer)
  FOREACH(FILE ${DGTAL_EXAMPLES_QGL_SRC}) 
    add_executable(${FILE} ${FILE})
    target_link_libraries ( ${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})   
    install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
  ENDFOREACH(FILE)
endif ( WITH_VISU3D_QGLVIEWER   )



add_executable(3dImageViewer 3dImageViewer specificClasses/Viewer3DImage  )
target_link_libraries ( 3dImageViewer   ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})   
