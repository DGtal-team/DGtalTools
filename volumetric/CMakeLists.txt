
SET(DGTAL_TOOLS_SRC
  criticalKernelsThinning3D
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
  volMask
)


FOREACH(FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE})
ENDFOREACH(FILE)


if (  WITH_VISU3D_QGLVIEWER )
  SET(DGTAL_EXAMPLES_QGL_SRC
    homotopicThinning3D )
  FOREACH(FILE ${DGTAL_EXAMPLES_QGL_SRC})
    DGtalTools_add_tool(${FILE})
  ENDFOREACH(FILE)
endif ( WITH_VISU3D_QGLVIEWER )
