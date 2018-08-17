# DGtalTools 0.9.4

- *converters*

   - mesh2vol: add option to add margin in the generated volume
     (to better extract the surfel boudary near domain limits).  
     (Bertrand Kerautret, [#322](https://github.com/DGtal-team/DGtalTools/pull/322))
   - vol2vox/vox2vol: tools to convert vol file to a MagicaVoxel VOX file and
   conversly. (David  Coeurjolly,
   [#314](https://github.com/DGtal-team/DGtalTools/pull/314))
   - volAddNoise moved to ```volumetric/```. (David  Coeurjolly,
   [#300](https://github.com/DGtal-team/DGtalTools/pull/300))
   - segfault fix in volBoundary2obj (David Coeurjolly,
   [#317](https://github.com/DGtal-team/DGtalTools/pull/317))
   - Fix the bad surfel display of volBoundary2obj (issue #320)
   (Bertrand Kerautret, [#321](https://github.com/DGtal-team/DGtalTools/pull/321))
   
- *volumetric*
   - new option to volAddNoise to extract the largest 6-connected
   component. (David  Coeurjolly,
   [#300](https://github.com/DGtal-team/DGtalTools/pull/300))
   - new option to 3dVolMarchingCubes to add some Kanungo noise to the
   input vol file. (David  Coeurjolly,
   [#300](https://github.com/DGtal-team/DGtalTools/pull/300))
   - Add thinning based of Critical Kernels using VoxelComplex.
   Based on [DGtal PR 1147](https://github.com/DGtal-team/DGtal/pull/1147)
   (Pablo Hernandez, [#311](https://github.com/DGtal-team/DGtalTools/pull/311))

- *visualisation*:
  - New tool for mesh voxelization from a mesh in input (.off)
    to a volumetric output (vol, pgm3d)
    (Monir Hadji, [#279](https://github.com/DGtal-team/DGtalTools/pull/279)
  - 2dCompImage : Computes and displays image comparisons (squared and absolute 
    differences) 
    (Bertrand Kerautret, [#313](https://github.com/DGtal-team/DGtalTools/pull/313))
  - Improve visualisation tools (vol2heightfield, vol2obj, vol2raw, vol2sdp, 
    vol2slice,volBoundary2obj,3dImageViewer, 3dVolViewer, sliceViewer, Viewer3DImage) 
    allowing to read longvol including rescaling. (Bertrand Kerautret, 
    [#296](https://github.com/DGtal-team/DGtalTools/pull/296))
  - Add an option to filter vector displayed in 3dSDPViewer.
   (Bertrand Kerautret, [#297](https://github.com/DGtal-team/DGtalTools/pull/297))

  - meshViewer: add an option to set the ambient light source.
    (Bertrand Kerautret, [#303](https://github.com/DGtal-team/DGtalTools/pull/303))
  - 3dSDPViewer: new option to display vector field as unit vectors.
    (Bertrand Kerautret, [#301](https://github.com/DGtal-team/DGtalTools/pull/304))

- *converters*:
  - sdp2vol: add the automatic set of the domain according to the
    bouding box of the set of points.    (Bertrand Kerautret,
    [#305](https://github.com/DGtal-team/DGtalTools/pull/305))



- *global*:
  - Fix travis Doxygen compilation for non Documention mode.
    (Bertrand Kerautret, [#314](https://github.com/DGtal-team/DGtalTools/pull/314))
  - Fix travis with boost installation which now use std package.
    (Bertrand Kerautret, [#310](https://github.com/DGtal-team/DGtalTools/pull/310))
  - Fix for the last QGLViewer version (2.7).
    (Bertrand Kerautret, [#308](https://github.com/DGtal-team/DGtalTools/pull/308))


# DGtalTools 0.9.3

- *global*:
   - Various fixes to enable the new Version3 (compressed) Vol/Longvol files.
     (David Coeurjolly, [#287](https://github.com/DGtal-team/DGtalTools/pull/287))
   - Fix Appveyor continuous integration with zlib installation and boost fix.
     (Bertrand Kerautret, [#289](https://github.com/DGtal-team/DGtalTools/pull/289))

- *imageProcessing*:
   - Creates imageProcessing directory. Add tools for doing image restoration
     and inpainting with Ambrosio-Tortorelli functional and discrete calculus.
     (Jacques-Olivier Lachaud, Marion Foare
     [#280](https://github.com/DGtal-team/DGtalTools/pull/280))

- *converters*:
   - fix tool itk2vol which was not able to read and convert int type image.
     (Bertrand Kerautret, [#276](https://github.com/DGtal-team/DGtalTools/pull/271))
   - add a CLOSURE export mode in volBoundary2obj. Default mode has been changed
     to "BDRY"
     (David Coeurjolly, [#281](https://github.com/DGtal-team/DGtalTools/pull/281))

- *visualisation*:
   - Add SnapShot option for meshViewer and 3dVolViewer
     (useful to get visualisation without interaction like for scripting and/or
     online demonstration). It also contains a new option to display a mesh in
     3DvolViewer.
     (Bertrand Kerautret, [#282](https://github.com/DGtal-team/DGtalTools/pull/282))
   - Add an option to display vector fields in displayContours
     (Bertrand Kerautret, [#290](https://github.com/DGtal-team/DGtalTools/pull/290))

- *estimators*:
    - 2dlocalEstimators: add an option to export the generated contour.
     (Bertrand Kerautret, [#285](https://github.com/DGtal-team/DGtalTools/pull/285))
    - tangentBC: add an option to read sdp points as input.
     (Bertrand Kerautret, [#285](https://github.com/DGtal-team/DGtalTools/pull/288))
    - volSurfaceRegularization: a tool to compute a regularized quadrangulation from
     from a digital surface.
     (Pierre Gueth, David Coeurjolly, [#306](https://github.com/DGtal-team/DGtalTools/pull/306))


# DGtalTools 0.9.2

- *global*:
  - fix wrong Khalimsky space initialization in Freeman2Img.
    (Roland Denis, [#271](https://github.com/DGtal-team/DGtalTools/pull/271))
  - doxygen documentation added for all tools. (David Coeurjolly, Bertrand Kerautret, [#258](https://github.com/DGtal-team/DGtalTools/pull/258))
  - fix uses of temporaries when ConstAlias is needed.
    (Roland Denis, [#253](https://github.com/DGtal-team/DGtalTools/pull/253))
  - renaming of the shapeGenerator folder to generators (David Coeurjolly, [#268](https://github.com/DGtal-team/DGtalTools/pull/268)))

- *visualisation*:
 - meshViewer: add a key to display mesh information about number of
    vertex/faces.
    (Bertrand Kerautret,
    [#273](https://github.com/DGtal-team/DGtalTools/pull/272))

  - 3dSDPViewer: fix the mesh display which was not given with their original
   colors. (Bertrand Kerautret,
   [#272](https://github.com/DGtal-team/DGtalTools/pull/272))

  - 3dSDPViewer: add the possibility to display a set of point by using
    different sphere sizes (specified in the input sdp file).
    (Bertrand Kerautret,
    [#252](https://github.com/DGtal-team/DGtalTools/pull/252))
  - sliceViewer: fix bug when imported image domain doesn't contain (0,0,0) point.
    (Roland Denis, [#256](https://github.com/DGtal-team/DGtalTools/pull/256))
  - 3dSDPViewer: add an option to display on screen the selected voxel.
    (Bertrand Kerautret,
    [#257](https://github.com/DGtal-team/DGtalTools/pull/257))


- *volumetric*:
  - fix reading options bug in volCComponentCounter and sdp2vol.
    (Bertrand Kerautret,
    [#254](https://github.com/DGtal-team/DGtalTools/pull/254))

# DGtalTools 0.9.1

- *converters*:
  - img2freeman: new option to sort the resulting contours by increasing size
    (B. Kerautret).

- *visualisation*:
  - meshViewer: new possibility to display a vector field from a simple sdp
    file (B. Kerautret).

  - 3dImplicitSurfaceExtractorByThickening: a tool to visualize 3d
    polynomial implicit surface defined as some f(x,y,z)=0. Its
    principle is to thickened the set {f=0} as {|f|<=e}, to extract an
    associated cubical complex, then to collapse it to capture the
    correct topology of {f=0}. Afterwards, the complex is projected
    onto f=0 with Newton's method. (J.-O. Lachaud)

  - 3dImplicitSurfaceExtractorBy4DExtension: a tool to visualize 3d
    polynomial implicit surface defined as some f(x,y,z)=0. Its
    principle is to extend f as a 4D function F(x,y,z,t)=0 (for
    instance F=f-|nabla f|t). This 4d hypersurface is easier to
    detect. It is transformed into 4D cubical complex, that is then
    collapsed to capture the correct topology of {f=0}. Afterwards,
    the complex is projected in 3D onto f=0 with Newton's
    method. (J.-O. Lachaud)

- *converters*:
  - homotopicThinning3D: the fixed points can be set from a file. (B. Kerautret)



# DGtalTools 0.9

- *global*:
  - Qt prog can now handle Qt5

- *visualisation*:
  - meshViewer: fix to display mesh with colored faces (.off with colors).
  - 3dCurvatureViewer/Noise: can now be used without QGLViewer if
    exporting data.
  - 3dCurvatureViewer/Noise:
    - can now be used without QGLViewer if exporting data
      (Jérémy Levallois).
    -  can now export II based normal vector field (Jérémy
      Levallois, David Coeurjolly).
  - 3dImageViewer: now display the current moving axis which was selected
	and display the slice numbers (can be disabled by key M).
  - sliceViewer: fix the bug when displaying vol with non 0 origin point.
  - itk2vol: convert any image of itk format (mhd, mha, ...) to vol
	(available with the itk option in DGtal).
  - sliceViewer: can now display 3d image with predefined color map (hue,
	gradient).

- *volumetric*:
  - volIntensityScale: a simple tool to apply a linear scale of the
    intensity given in a volumetric file.

- *converters*:
  - vol2heightfield: a new tool to transform volumetric file into 2D
    heightmap.
  - heightfield2vol: a new tool to transform 2D heightmap into volumetric
    file.
  - imgAddNoise: a new tool to add noise (Kanungo's) to a binary 2D object
  - volAddNoise: a new tool to add noise (Kanungo's) to a binary 3D object
  - heightfield2shading: a new tool to render a 2D heightfield image into
    a shading one.
  - mesh2heightfield:  new tool to convert a mesh file into a 2D heightmap
    (from a normal direction N and from a starting point P).
  - freeman2img: (extended from previous freeman2pgm) fix options issues.


- *estimators*:
  - 3dCurveTangentEstimator: a simple tool to estimate and visualize the tangent to a set of points approaching a 3D curve.
    Two estimators are implemented, one based on digital Voronoi Covariance Measure, the other based on 3D lambda-MST.

# DGtalTools 0.8

- *visualisation*:
  - 3dCompSurfelData: a tool to compare generic surfel data informations given from two data files.
  - 3dCurvatureViewer: can now display curvature on multiple connected components and can apply image re sampling for anisotropic grid.
  - 3dDisplaySurfelData: display surfel data from SDP file with color attributes given as scalar interpreted as color.
  - 3dHeightMapViewer: display a 2D image as heightmap by using QGLviewer.
  - 3dSDPViewer: basic display of a sequence of 3d points (as voxel or sphere) and vectors by using QGLviewer.
  - 3dVolBoundaryViewer: a simple viewer of the boundary of an object defined by thresholding a vol file.
  - sliceViewer: a new 2D and 3D slice viewer from 3D volumic files ( pgm3d, vol, longvol, and DICOM with ITK).

- *converters*:
  - freeman2pgm: transform one or several freeman chains into a pgm file by filling their interior areas and renamed into freeman2img.
  - pgm2freeman: renamed to pgm2img.
  - vol2off: tool removed, see 3dVolMarchingCubes for the same tool.
  - volBoundary2obj: a simple tool to export the boundary of a an	object in a volumetric file to OBJ

- *estimators*:
  - curvatureScaleSpaceBCC: a tool to display the curvature scale space of a given contour with the Binomial Convolver Curvature Estimator.
  - eulerCharacteristic: bruteforce tool to extract (volumetric) Euler characteristic from volumetric binary object.
  - generic3dNormalEstimators: Computes a normal vector field over a digitized 3D implicit surface for several estimators (II|VCM|Trivial|True).

- *volumetric*:
  - 3dVolMarchingCubes: speed-up by factor 10 simply by replacing the set predicate used in computations.
  - volCrop: crop an 3D vol image from to points.
  - volReSample: apply a basic  re sampling of a 3D volumetric image (.vol, .longvol, .pgm3d)  with a given grid size.
  - volSegment: Segment volumetric file from a simple threshold which can be set automatically from the Otsu's variance based estimation.
  - volTrValues: a basic tool to transform the voxel values from an input/output set.


# DGtalTools 0.7

- *converters*:
  - convertVol: a simple generic volume image converters (can process actually pgm3d, vol, longvol, raw (for writing)).
  - vol2sdp: a simple tools to extract digital points from 3d vol files.
  - vol2off: extract dual surface of a digital object (equiv. Marching Cubes)
  - vol2obj: convert a volume file into OBJ format (all voxels belonging to threshold interval)
  - vol2slice: tool to extract all slices from 3d volumic images.
  - slice2vol: tool to merge slices into one 3d volumic file.
  - sdp2vol: a simple tool to create a 3d vol image from 3d digital points.
  - longvol2vol: convert longvol to vol file using different conversion policies.
  - dicom2vol: convert dicom images into 3d volumic file (need itk option in DGtal).
  - pgm2freeman: add new possibility to set automatically a threshold from the otsu algorithm.
  - HDF52vol: convert HDF5 to vol file format.
  - raw2HDF5: convert raw image to HDF5.

- *volumetric*:
  - homotopicThinning3D exploits now the GenericReader class and is no more limited to vol format.
  - volFlip: tool to flip all volume slice images according a given dimension.
  - volImageMetrics: apply basic measures from two volumetric images:  RMSE and PSNR.
  - volShapeMetrics: Apply shape measures for comparing two volumetric images A and B (shape defined from thresholds):
    - Measures from voxel partition (true/false+-, precision recall, f-measure)
    - Measures bases on euclidean distance between the two Shape A and B.

- *estimators*:
  - 2dLocalEstimators: Improvement of 2dLocalEstimators + possibility to compare with noised data.
  - 3dLocalEstimators: Adding possibility to compare curvature (mean, gaussian and principal curvatures)
     with Integral Invariant and Monge via Jet Fitting + possibility to compare with noised data.

- *volumetric*:
  - volTools directory moved into volumetric.

- *visualisation*:
  - 3dCurveViewer: A tool for visualizing the tangential cover of 3d curves.
  - 3dVolViewer: new option to limit the number of displayed voxels (can open dicom format if WITH_ITK is set to true).
  - 3dImageViewer: new tool to display slice image with interactive translatations or rotations (can open dicom format  if WITH_ITK is set to true).
  - patternTriangulation: a new tool that draws with Board2D the convex hull, the closest-point Delaunay triangulation or the farthest-point Delaunay triangulation of a pattern.
  - 3dCurvatureViewer: Now allow to draw principal curvature directions on objets.
  - 3dCurvatureViewerNoise: Same as 3dCurvatureViewer, but allow to add some noise to objects.


- *distanceTransform*:
  - LUTBasedNSDistanceTransform: Compute the 2D translated neighborhood-sequence distance transform of a binary image.
  - CumulativeSequenceTest and RationalBeattySequenceTest: tests from LUTBasedNSDistanceTransform.

# DGtalTools 0.6

- *estimators*:
  - 2dLocalEstimators: program to compare local curvature/tangent estimators on implicit shapes:
    - Maximal DSS based estimators
    - Maximal DCA based estimators
    - Binomial convolver based estimators
    - Integral Invariants based estimators
  -3dLocalEstimators: program to compare  3D local curvature (mean or gaussian) estimators on 3D implicit shapes.

- *visualisation*:
  - 3dCurvatureViewer: computes and displays mean or gaussian curvature of vol binary shapes.
   - Various updates for 0.6 DGtal compatibility.




# DGtalTools 0.1

- *converters*: utilities to convert various simple file formats:
  - freeman2sdp: convert freeman chain towards a Sequence of Discrete Points.
  - pgm2freeman: to extract a freeman chain contour from a grayscale image.
  - raw2vol and vol2raw: transform 3D volumes files from (resp. to) raw to vol.
  - ofs2off: convert OFS mesh format towards a OFF variant.

- *estimators*:
  - lengthEstimator: program to generate multigrid analysis of length estimators.
  - tangentBC: tangent estimator using the Binomial convolver.
  - curvatureBC: curvature estimator using the Binomial convolver.
  - curvatureMCMS: curvature estimator using the maximal segments cover  (to be updated for current DGtal version).
  - estimatorComparator: program to perform comparison of local quantity estimators (to be updated for current DGtal version).
  - vol2normalField: compute the normal vector field of a given vol file .

- *shapeGenerator*:
  - shapeGenerator: generate multigrid shape
  - contourGenerator: generate multigrid shape contours


- *visualization*:
  - 3dVolViewer: volume file (.vol and .pgm3d) viewer with QGLViewer.
  - displayContours: display discrete contours from various format (.fc (freemanchain), .sdp).
  - meshViewer: display 3D mesh from OFS or OFF format.

- *volumetric*:
  - 3dVolMarchingCubes: marching cubes form a Vol file
  - homotopicThinning3D: ultimate skeleton from vol file

- *volumetric:
  - volAddBorder: add a 1 voxel boundary with value 0 to a vol file.
  - volCComponentCounter: a simple program to count the number of connected components in a 3D image.
  - volSubSample: sub sample a vol file (division by 2 in each direction).
