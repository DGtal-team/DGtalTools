SET(DGTAL_TOOLS_SRC
	displayContours
        patternTriangulation
        2dCompImage
)

foreach (FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE})
ENDFOREACH(FILE)

if (  DGTAL_WITH_POLYSCOPE_VIEWER )
    SET(DGTAL_EXAMPLES_POLYSCOPE_SRC
    3dVolViewer
    3dVolBoundaryViewer
    3dImageViewer
    3dCurveViewer
    3dCurvatureViewer
    3dCurvatureViewerNoise
    3dSDPViewer
    meshViewer
    3dHeightMapViewer
    3dDisplaySurfelData
    3dCompSurfelData
    3dImplicitSurfaceExtractorBy4DExtension
    3dImplicitSurfaceExtractorByThickening )

  FOREACH(FILE ${DGTAL_EXAMPLES_POLYSCOPE_SRC})
    DGtalTools_add_tool(${FILE})
  endforeach (FILE)


  include_directories(${CMAKE_CURRENT_BINARY_DIR})
endif (  DGTAL_WITH_POLYSCOPE_VIEWER )


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



