
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
  volFillInterior
  volInfo
  volMask
)


FOREACH(FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE})
ENDFOREACH(FILE)


if ( DGTAL_WITH_POLYSCOPE_VIEWER )
  SET(DGTAL_EXAMPLES_POLYSCOPE_SRC
    criticalKernelsThinning3D
    homotopicThinning3D )
  FOREACH(FILE ${DGTAL_EXAMPLES_POLYSCOPE_SRC})
    DGtalTools_add_tool(${FILE})
  ENDFOREACH(FILE)
endif ( DGTAL_WITH_POLYSCOPE_VIEWER )
