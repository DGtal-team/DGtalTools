
SET(DGTAL_TOOLS_SRC
  imgAddNoise
  img2freeman
  freeman2img
  freeman2sdp
  vol2raw
  raw2vol
  ofs2off
  convertVol
  vol2sdp
  slice2vol
  vol2slice
  vol2obj
  longvol2vol
  sdp2vol
  volBoundary2obj
  vol2heightfield
  heightfield2vol
  heightfield2shading
  mesh2heightfield
  vox2vol
  vol2vox
  mesh2vol)

if( WITH_HDF5 )
SET(DGTAL_TOOLS_SRC ${DGTAL_TOOLS_SRC}
    raw2HDF5
    HDF52vol)
endif( WITH_HDF5 )

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  DGtalTools_add_tool(${FILE})
ENDFOREACH(FILE)


if (  WITH_ITK )
  SET(DGTAL_TOOLS_ITK
    dicom2vol
    itk2vol)
  FOREACH(FILE ${DGTAL_TOOLS_ITK})
    DGtalTools_add_tool(${FILE})
  ENDFOREACH(FILE)
ENDIF( WITH_ITK )
