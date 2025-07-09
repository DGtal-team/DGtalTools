SET(DGTAL_TOOLS_SRC
  2dLocalEstimators
  lengthEstimators
  vol2normalField
  curvatureBC
  curvatureMCMS
  tangentBC
  statisticsEstimators
  generic3dNormalEstimators
  curvatureScaleSpaceBCC
  eulerCharacteristic
)

if (DGTAL_WITH_EIGEN)
    DGtalTools_add_tool(volSurfaceRegularization volSurfaceRegularization-details/surface_approx.cpp)
endif()

if ( DGTAL_WITH_POLYSCOPE_VIEWER )
  SET(VISU3D_TESTS_SRC
      3dCurveTangentEstimator
  )
  FOREACH(FILE ${VISU3D_TESTS_SRC})
     DGtalTools_add_tool(${FILE})
  ENDFOREACH(FILE)
endif ( DGTAL_WITH_POLYSCOPE_VIEWER )

if (  DGTAL_WITH_CGAL )
  SET(CGAL_TESTS_SRC
    3dLocalEstimators
     )
  FOREACH(FILE ${CGAL_TESTS_SRC})
     DGtalTools_add_tool(${FILE})
  ENDFOREACH(FILE)
endif ( DGTAL_WITH_CGAL )

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE})
  if(MSVC)
    set_source_files_properties( ${FILE} PROPERTIES COMPILE_FLAGS "/bigobj" )
  endif(MSVC)
ENDFOREACH(FILE)
