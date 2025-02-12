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

if (WITH_EIGEN)
    DGtalTools_add_tool(volSurfaceRegularization volSurfaceRegularization-details/surface_approx.cpp)
endif()

if ( WITH_VISU3D_QGLVIEWER )
  SET(VISU3D_TESTS_SRC
      3dCurveTangentEstimator
  )
  FOREACH(FILE ${VISU3D_TESTS_SRC})
     DGtalTools_add_tool(${FILE})
  ENDFOREACH(FILE)
endif ( WITH_VISU3D_QGLVIEWER )

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
