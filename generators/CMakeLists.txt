SET(DGTAL_TOOLS_SRC
  shapeGenerator
  contourGenerator
  3dParametricCurveDigitizer
)

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE})
ENDFOREACH(FILE)

