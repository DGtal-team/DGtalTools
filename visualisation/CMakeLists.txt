SET(DGTAL_TOOLS_SRC
	displayContours
        patternTriangulation
        3dCurvatureViewer
        3dCurvatureViewerNoise
        2dCompImage
)

foreach (FILE ${DGTAL_TOOLS_SRC})
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
    3dCurveViewer
    3dSDPViewer
    meshViewer
    3dHeightMapViewer
    3dDisplaySurfelData
    3dCompSurfelData
    3dImplicitSurfaceExtractorBy4DExtension
    3dImplicitSurfaceExtractorByThickening )

  FOREACH(FILE ${DGTAL_EXAMPLES_QGL_SRC})
    add_executable(${FILE} ${FILE})
    target_link_libraries ( ${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})

    install(TARGETS ${FILE} RUNTIME DESTINATION bin
                            LIBRARY DESTINATION lib
                            ARCHIVE DESTINATION lib)
  endforeach (FILE)

  add_executable(3dImageViewer 3dImageViewer specificClasses/Viewer3DImage)
  target_link_libraries(3dImageViewer ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
  install(TARGETS 3dImageViewer RUNTIME DESTINATION bin
                                LIBRARY DESTINATION lib
                                ARCHIVE DESTINATION lib)

  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  if (WITH_QT4)
    QT4_WRAP_UI(ui_sliceViewer_FORMS_HEADERS sliceViewer.ui)
    QT4_WRAP_CPP(sliceViewer_HEADERS_MOC sliceViewer.h OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)

    add_executable(sliceViewer sliceViewer.cpp ${sliceViewer_HEADERS_MOC}
                                               ${ui_sliceViewer_FORMS_HEADERS})
    target_link_libraries (sliceViewer ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})

  elseif (WITH_QT5)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    qt5_wrap_ui(ui_sliceViewer_FORMS_HEADERS sliceViewer.ui)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5OpenGL_EXECUTABLE_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Xml_EXECUTABLE_COMPILE_FLAGS}")
    add_executable(sliceViewer sliceViewer.cpp ${ui_sliceViewer_FORMS_HEADERS})
    target_link_libraries(sliceViewer ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})
    qt5_use_modules(sliceViewer Widgets OpenGL Xml)
  endif (WITH_QT4)

  install(TARGETS sliceViewer RUNTIME DESTINATION bin
                              LIBRARY DESTINATION lib
                              ARCHIVE DESTINATION lib)
endif (WITH_VISU3D_QGLVIEWER)



