
SET(DGTAL_TOOLS_SRC
  volAddNoise
  3dVolMarchingCubes
  volAddBorder
  volCComponentCounter
  volSubSample
  volCrop
  volImageMetrics
  volShapeMetrics
  volFlip
  volSegment
  volReSample
  volTrValues
  volIntensityScale
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
    homotopicThinning3D
    criticalKernelsThinning3D )
  FOREACH(FILE ${DGTAL_EXAMPLES_QGL_SRC})
    add_executable(${FILE} ${FILE})
    target_link_libraries ( ${FILE} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
    install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
  ENDFOREACH(FILE)
endif ( WITH_VISU3D_QGLVIEWER   )
