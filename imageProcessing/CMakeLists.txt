
# These tools require Eigen (because DEC requires Eigen)
if(DGTAL_WITH_EIGEN)
  SET(DGTAL_TOOLS_SRC
    at-u0-v1
    at-u2-v0
    atv-u2-v0
    )
  FOREACH(FILE ${DGTAL_TOOLS_SRC})
    DGtalTools_add_tool(${FILE})
  ENDFOREACH(FILE)
endif(DGTAL_WITH_EIGEN)
