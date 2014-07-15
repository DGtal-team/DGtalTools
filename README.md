DGtalTools
==========

DGtalTools is a separate github project containing tools constructed
using DGtal library. The main goal of this part is to gather simple
and useful tools exploiting the structures and algorithms defined in
DGtal. The resulting tools could be useful to:

- Share and apply DGtal algorithms to various data from different domains.

- Construct demonstration tools like online demonstrations (as for instance the one of the Image Processing Online (http://www.ipol.im)

- Simplify comparisons of different algorithms with an single framework.

- Provide useful tools of digital image related algorithms (extraction
  of connected components, digital contour/surface extraction, simple
  visualization tools ... etc).



The source code of the tools can also be used to non DGtal familiar
user to show how to include the DGtal library framework directly in their
own source code (in complement of DGtal tutorial http://libdgtal.org/doc/stable/tutorials.html).


Organisation
============

Actually the DGTal project is organized as follows:

 - [Converters](#converters)

 - [DistanceTransform] (#DistanceTransform)

 - [Estimators](#estimators)

 - [ShapeGenerator](#shapeGenerator)

 - [Visualization](#visualizatoin)

 - [Volumetric](#volumetric)


Converters
----------

Utilities to convert various simple file formats:

  - convertVol: a simple generic volume image converters (can process actually pgm3d, vol, longvol, raw (for writing)).
  - dicom2vol: convert dicom images into 3d volumic file.
  - freeman2sdp: convert freeman chain towards a Sequence of Discrete Points.
  - HDF52vol: convert HDF5 to vol file format. 
  - longvol2vol: convert longvol to vol file using different conversion policies.  
  - ofs2off: convert OFS mesh format towards a OFF variant.
  - pgm2freeman: to extract a freeman chain contour from a grayscale image.
  - freeman2pgm: transform one or several freeman chains into a pgm file by filling their interior areas.
  - raw2vol and vol2raw: transform 3D volumes files from (resp. to) raw to vol.
  - raw2HDF5: convert raw image to HDF5.
  - slice2vol: tool to merge slices into one 3d volumic file.
  - sdp2vol: a simple tool to create a 3d vol image from 3d digital points.
  - vol2sdp: a simple tools to extract digital points from 3d vol files.
  - vol2off: extract dual surface of a digital object (equiv. Marching Cubes)
  - vol2obj: convert a volume file into OBJ format (all voxels belonging to threshold interval)   
  - vol2slice: tool to extract all slices from 3d volumic images.
  
DistanceTransform
------------------

  - LUTBasedNSDistanceTransform: Compute the 2D translated neighborhood-sequence distance transform of a binary image.
  - CumulativeSequenceTest and RationalBeattySequenceTest: tests from LUTBasedNSDistanceTransform.



Estimators
----------
  - generic3dNormalEstimators: Computes a normal vector field over a digitized 3D implicit surface for several estimators (II|VCM|Trivial|True).
  - 2dLocalEstimators: program to compare local curvature/tangent estimators on implicit shapes
    - Maximal DSS based estimators
    - Maximal DCA based estimators
    - Binomial convolver based estimators
    - Integral Invariants based estimators
  - 3dLocalEstimators: program to compare  3D local curvature (mean or gaussian) estimators on 3D implicit shapes.
  - lengthEstimator: program to generate multigrid analysis of length estimators.
  - tangentBC: tangent estimator using the Binomial convolver.
  - curvatureBC: curvature estimator using the Binomial convolver.
  - curvatureMCMS: curvature estimator using the maximal segments cover  (to be updated for current DGtal version).
  - curvatureScaleSpaceBCC: a tool to display the curvature scale space of a given contour with the Binomial Convolver Curvature Estimator
  - estimatorComparator: program to perform comparison of local quantity estimators (to be updated for current DGtal version).
  - vol2normalField: compute the normal vector field of a given vol file .
<table>
<tr>
<td colspan="2"><img height=130 src="https://cloud.githubusercontent.com/assets/772865/2646108/f515b0a2-bf39-11e3-96f8-c7606173f43b.png"></td>
</tr>
<tr>
<td colspan="2">Illustration of curvatureScaleSpaceBCC </td>
</tr>
<tr>
<td align=center ><img height=130 src="https://cloud.githubusercontent.com/assets/793707/2996392/d3ee9e58-dced-11e3-98a0-72233927aaf6.jpg"> </td>
<td align=center ><img height=130 src="https://cloud.githubusercontent.com/assets/772865/3311642/03f29044-f6c7-11e3-8981-120369b3e8bd.png"> </td>
</tr>
<tr>
<td colspan="2"> Illustration of generic3dNormalEstimators on VCM estimator applied on smooth and noisy shapes.</td>
</tr>

</table>

ShapeGenerator
--------------
  - shapeGenerator: generate multigrid shape
  - contourGenerator: generate multigrid shape contours
 <center>
<table>
<tr>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684690/eff46c16-da02-11e2-861e-ddc366b247e8.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684694/39a9cc2a-da03-11e2-9f49-3aff0e886c35.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684695/42b657ca-da03-11e2-985e-e468084b5c01.png"></td>
</tr>
<tr>
<td> grid size = 1</td> <td> grid size= 0.1</td> <td> grid size = 0.01</td>
</tr>
<tr>
<td colspan=3 > Illustration of the shapeGenerator tools </td>
</tr>
</table>
</center>

Visualization
-------------
  - 3dCurvatureViewer: permits to compute and visualize mean or gaussian curvature of binary shapes.
  - 3dCurvatureViewerNoise: Same as 3dCurvatureViewer, but allow to add some noise to objects.
  - 3dCurveViewer: displays 3D curves with tangential cover and projections onto bounding box.
  - 3dImageViewer: tools to display 3d slice images (.vol, .pgm3d and  dicom with ITK) with QGLViewer.
  - 3dSDPViewer: basic display of a sequence of 3d points (as voxel or sphere) and vectors by using QGLviewer.
  - 3dVolViewer: volume file (.vol, .pgm3d and dicom with ITK ) viewer with QGLViewer.
  - displayContours: display discrete contours from various format (.fc (freemanchain), .sdp).
  - meshViewer: display 3D mesh from OFS or OFF format.
  - patternTriangulation: a new tool that draws with Board2D the convex hull, the closest-point Delaunay triangulation or the farthest-point Delaunay triangulation of a pattern.


Here are some illustrations of such a tools:
<center>
<table>
<tr>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684607/450c064e-da00-11e2-8830-76eb90a5efd7.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/685853/d96a5252-da44-11e2-9872-7f0160be8f5d.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684569/59a2f6fa-d9fe-11e2-84ba-a48842f4aafb.png" ></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684590/778bea9a-d9ff-11e2-8e04-6e3e8a39ae3c.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684598/c3adcf4c-d9ff-11e2-8c3f-e67c8abd0c76.png"></td>
</tr>
<tr>
<td>3dCurvatureViewer</td>
<td>3dCurveViewer </td>
<td>3dImageViewer</td>
<td>3dVolViewer</td>
<td>displayContours</td>
</tr>
<tr>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684622/d698405a-da00-11e2-8aa0-19212a58ce23.png"></td>
<td><img height=130 src="https://cloud.githubusercontent.com/assets/772865/2720141/6c42a0e0-c56b-11e3-8328-a6d88242f21e.png"> </td>
</tr>
<tr>
<td>meshViewer</td>
<td>3dSDPViewer</td>

</tr>



</table>
</center>
Volumetric
----------
  - 3dVolMarchingCubes: marching cubes form a Vol file
  - homotopicThinning3D: ultimate skeleton from vol file
  - volAddBorder: add a 1 voxel boundary with value 0 to a vol file.
  - volCComponentCounter: a simple program to count the number of connected components in a 3D image.
  - volFlip: tool to flip all volume slice images according a given dimension.
  - volSubSample: sub sample a vol file (division by 2 in each direction).
  - volImageMetrics: apply basic statistics on comparaison between two volumetric images (shape defined from thresholds): computes true/false -+, precision, recall f-mean RMSE, PSNR.
  - volShapeMetrics: apply euclidean distance comparisons between two shapes  (shape defined from thresholds).
  - volReSample: apply a basic  re sampling of a 3D volumetric image (.vol, .longvol, .pgm3d)  with a given grid size.


How to build the tools
======================
  - use cmake tool to generate a build script (MakeFile, VS project,..) from the CMakeLists.txt
  - DGtal must be installed in your system. Concerning DGtal dependencies (boost, Qt,...), all the dependencies used to compile your DGtal library must be present to build the DGtalTools.
  
  
