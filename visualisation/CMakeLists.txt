SET(DGTAL_TOOLS_SRC
	displayContours
        patternTriangulation
        3dCurvatureViewer
        3dCurvatureViewerNoise
        2dCompImage
)

foreach (FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE})
ENDFOREACH(FILE)

if (  WITH_VISU3D_QGLVIEWER )
    set_property(SOURCE ui_sliceViewer.h PROPERTY SKIP_AUTOMOC ON)
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
    DGtalTools_add_tool(${FILE})
  endforeach (FILE)

  DGtalTools_add_tool(3dImageViewer specificClasses/Viewer3DImage.cpp)

  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  if (WITH_QT4)
    QT4_WRAP_UI(ui_sliceViewer_FORMS_HEADERS sliceViewer.ui)
    QT4_WRAP_CPP(sliceViewer_HEADERS_MOC sliceViewer.h OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
      DGtalTools_add_tool(sliceViewer ${sliceViewer_HEADERS_MOC} ${ui_sliceViewer_FORMS_HEADERS})      
  elseif (WITH_QT5)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    qt5_wrap_ui(ui_sliceViewer_FORMS_HEADERS sliceViewer.ui)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5OpenGL_EXECUTABLE_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Xml_EXECUTABLE_COMPILE_FLAGS}")
    DGtalTools_add_tool(sliceViewer ${ui_sliceViewer_FORMS_HEADERS})
  endif (WITH_QT4)
endif (WITH_VISU3D_QGLVIEWER)



if (  NOT(DISABLE_POLYSCOPE) )
    SET(DGTALTOOLS_POLYSCOPE_SRC
    volscope)

  FOREACH(tool_file ${DGTALTOOLS_POLYSCOPE_SRC})
    add_executable(${tool_file} ${tool_file}.cpp)
    target_link_libraries (${tool_file} ${DGTAL_LIBRARIES} ${DGtalLibDependencies} ${DGtalToolsLibDependencies} polyscope)
    install(TARGETS ${tool_file}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
  endforeach (tool_file)
endif ()



