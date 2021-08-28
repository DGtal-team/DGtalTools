
SET(DGTAL_TOOLS_SRC
 rigidTransform2D
 rigidTransform3D
)


FOREACH(FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE} ${FILE}.cpp)
ENDFOREACH(FILE)

