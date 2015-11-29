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
own source code (in complement of DGtal tutorial http://liris.cnrs.fr/dgtal/doc/nightly/packageTutorials.html).


* Release 0.9 DOI [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.31882.svg)](http://dx.doi.org/10.5281/zenodo.31882)
* Release 0.8 DOI [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11587.svg)](http://dx.doi.org/10.5281/zenodo.11587)
* [![Build Status](https://api.travis-ci.org/DGtal-team/DGtalTools.svg?branch=master)](https://travis-ci.org/DGtal-team/DGtalTools)


Organisation
============

Actually the DGTal project is organized as follows:

 - [Converters](#converters)

 - [DistanceTransform](#distancetransform)

 - [Estimators](#estimators)

 - [ShapeGenerator](#shapegenerator)

 - [Visualization](#visualization)

 - [Volumetric](#volumetric)


Converters
----------

Utilities to convert various simple file formats:

  - convertVol: a simple generic volume image converters (can process actually pgm3d, vol, longvol, raw (for writing)).
  - dicom2vol: convert dicom images into 3d volumic file.
  - freeman2img: transform one or several freeman chains into a pgm file by filling their interior areas.
  - freeman2sdp: convert freeman chain towards a Sequence of Discrete Points.
  - HDF52vol: convert HDF5 to vol file format. 
  - heightfield2shading: Render a 2D heightfield image into a shading image from various reflectance models (lambertian, specular, custom reflectance map).
  - heightfield2vol: a new tool to transform 2D heightmap into volumetric file.
  - img2freeman: to extract a freeman chain contour from a grayscale image.
  - imgAddNoise: a new tool to add noise (Kanungo's) to a binary 2D object.
  - itk2vol: convert any image of itk format (mhd, mha, ...) to vol (available with the itk option in DGtal).
  - longvol2vol: convert longvol to vol file using different conversion policies.  
  - mesh2heightfield: a tool to convert a mesh file into a 2D heightmap (from a normal direction N and from a starting point P).
  - ofs2off: convert OFS mesh format towards a OFF variant.
  - raw2HDF5: convert raw image to HDF5.
  - raw2vol and vol2raw: transform 3D volumes files from (resp. to) raw to vol.
  - sdp2vol: a simple tool to create a 3d vol image from 3d digital points.
  - slice2vol: tool to merge slices into one 3d volumic file.
  - volAddNoise: a new tool to add noise (Kanungo's) to a binary 3D object.
  - vol2obj: convert a volume file into OBJ format (all voxels belonging to threshold interval)   
  - vol2raw: convert a vol to a 8-bit raw file.
  - vol2sdp: a simple tools to extract digital points from 3d vol files.
  - vol2slice: tool to extract all slices from 3d volumic images.
  - volBoundary2obj: export the boundary of a volume file to OBJ format.
  - vol2heightfield: a new tool to transform volumetric file into 2D heightmap.

<center>
<table>
<tr>
<td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/6706730/9bac9720-cd60-11e4-9819-81e536b21e97.gif"> </td>
<td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/6706785/0ca0071e-cd61-11e4-9304-c6e168b1c6b2.png"> </td>
<td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/6914017/ccf433e4-d786-11e4-997b-f513f07f56f3.gif"> </td>
<td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/6914104/6d25be78-d787-11e4-9433-2a834fc4a0af.png"> </td>
</tr>
<tr>
<td>vol2heightmap</td>
<td>heightmap2vol</td>
<td>imgAddNoise</td> 
<td>volAddNoise</td>
</tr>
<tr>
<td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/7106643/c09cb8d4-e148-11e4-8653-2d5bac5dc3c5.png"> 
<td align=center> <img height=150 src= "https://cloud.githubusercontent.com/assets/772865/10265156/64aaad64-6a24-11e5-9773-628c661fa76c.png">
</td>
</tr>
<tr>
<td>mesh2heightfield</td>
<td>heightfield2shading</td>
</tr>

<tr>

</table>

</center>




DistanceTransform
-----------------
  - LUTBasedNSDistanceTransform: Compute the 2D translated neighborhood-sequence distance transform of a binary image.
  - CumulativeSequenceTest and RationalBeattySequenceTest: tests from LUTBasedNSDistanceTransform.
  


Estimators
----------

  - 2dLocalEstimators: program to compare local curvature/tangent estimators on implicit shapes
    - Maximal DSS based estimators
    - Maximal DCA based estimators
    - Binomial convolver based estimators
    - Integral Invariants based estimators
  - 3dLocalEstimators: program to compare  3D local curvature (mean or gaussian) estimators on 3D implicit shapes.
  - curvatureBC: curvature estimator using the Binomial convolver.
  - curvatureMCMS: curvature estimator using the maximal segments cover  (to be updated for current DGtal version).
  - curvatureScaleSpaceBCC: a tool to display the curvature scale space of a given contour with the Binomial Convolver Curvature Estimator
  - eulerCharacteristic: bruteforce tool to extract (volumetric) Euler characteristic from volumetric binary object.
  - generic3dNormalEstimators: computes a normal vector field over a digitized 3D implicit surface for several estimators (II|VCM|Trivial|True).
  - lengthEstimators: generate multigrid length estimations of paramteric shapes using DGtal library.
  - statisticsEstimators: compute satistics (L1, L2, Loo) from results of two extimators.
  - tangentBC: tangent estimator using the Binomial convolver.
  - vol2normalField: compute the normal vector field of a given vol file.

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
  - contourGenerator: generate multigrid shape contours.
  - shapeGenerator: generate multigrid shape.
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
  - 3dCompSurfelData: it computes generic scalar surfel data comparisons with display with surfel association.
  - 3dCurvatureViewer: permits to compute and visualize mean or gaussian curvature of binary shapes.
  - 3dCurvatureViewerNoise: Same as 3dCurvatureViewer, but allow to add some noise to objects.
  - 3dCurveViewer: displays 3D curves with tangential cover and projections onto bounding box.
  - 3dDisplaySurfelData: displays surfel data from SDP file with color attributes given as scalar interpreted as color.
  - 3dHeightMapViewer: display a 2D image as heightmap by using QGLviewer.
  - 3dImageViewer: tools to display 3d slice images (.vol, .pgm3d and  dicom with ITK) with QGLViewer.
  - 3dSDPViewer: basic display of a sequence of 3d points (as voxel or sphere) and vectors by using QGLviewer.
  - 3dVolBoundaryViewer: Display the boundary of a volume file by using QGLviewer. 
  - 3dVolViewer: volume file (.vol, .pgm3d and dicom with ITK ) viewer with QGLViewer.
  - displayContours: display discrete contours from various format (.fc (freemanchain), .sdp).
  - meshViewer: display 3D mesh from OFS or OFF format.
  - patternTriangulation: a new tool that draws with Board2D the convex hull, the closest-point Delaunay triangulation or the farthest-point Delaunay triangulation of a pattern.
  - sliceViewer: a new 2D and 3D slice viewer from 3D volumic files ( pgm3d, vol, longvol, and DICOM with ITK).

Here are some illustrations of such a tools:
<center>
<table>
<tr>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684607/450c064e-da00-11e2-8830-76eb90a5efd7.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/685853/d96a5252-da44-11e2-9872-7f0160be8f5d.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684569/59a2f6fa-d9fe-11e2-84ba-a48842f4aafb.png" ></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684590/778bea9a-d9ff-11e2-8e04-6e3e8a39ae3c.png"></td>
</tr>
<tr>
<td>3dCurvatureViewer</td>
<td>3dCurveViewer </td>
<td>3dImageViewer</td>
<td>3dVolViewer</td>
</tr>

<tr>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684598/c3adcf4c-d9ff-11e2-8c3f-e67c8abd0c76.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684622/d698405a-da00-11e2-8aa0-19212a58ce23.png"></td>
<td><img width=300 src="https://cloud.githubusercontent.com/assets/772865/2720141/6c42a0e0-c56b-11e3-8328-a6d88242f21e.png"> </td>
<td><img  width=300 src="https://cloud.githubusercontent.com/assets/772865/10269635/8678ca02-6add-11e5-9d83-fbcf3608612f.png"></td>
</tr>
<tr>
<td>displayContours</td>
<td>meshViewer</td>
<td>3dSDPViewer</td>
<td>sliceViewer</td>
</tr>

<tr>
<td><img  width=300 src="https://cloud.githubusercontent.com/assets/772865/3486505/6edb2144-043e-11e4-81c4-2c20f272a119.png" </td>
<td><img width=300 src="https://cloud.githubusercontent.com/assets/772865/4118303/5f53f61c-329f-11e4-9629-23c53afd9eff.png" </td>
</tr>
<tr>
<td>3dHeightMapViewer</td>
<td>3dCompSurfelData</td>
</tr>



</table>
</center>
Volumetric
----------
  - 3dVolMarchingCubes: marching cubes form a Vol file.
  - homotopicThinning3D: ultimate skeleton from vol file.
  - volAddBorder: add a 1 voxel boundary with value 0 to a vol file.
  - volCComponentCounter: a simple program to count the number of connected components in a 3D image.
  - volCrop: Crop an 3D vol image from to points.
  - volFlip: tool to flip all volume slice images according a given dimension.
  - volImageMetrics: apply basic statistics on comparaison between two volumetric images (shape defined from thresholds): computes true/false -+, precision, recall f-mean RMSE, PSNR. 
  - volReSample: apply a basic  re sampling of a 3D volumetric image (.vol, .longvol, .pgm3d)  with a given grid size. 
  - volSegment: Segment volumetric file from a simple threshold which can be set automatically from the otsu estimation.
  - volShapeMetrics: apply euclidean distance comparisons between two shapes  (shape defined from thresholds).
  - volSubSample: sub sample a vol file (division by 2 in each direction).
  - volTrValues: apply basic vol image transform from the input values to output values.

<center>
<table>
<tr>
<td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/4140989/2b115d34-33aa-11e4-9141-749d02335a4d.png"> </td>
<td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/4140972/f13c7fe4-33a9-11e4-8114-74e4925f0628.png"> </td>
</tr>
<tr>
<td>homotopicThinning3D</td>
<td>volSubSample</td>
</table>
</center>



How to build the tools
======================
  - use cmake tool to generate a build script (MakeFile, VS project,..) from the CMakeLists.txt
  - DGtal must be installed in your system. Concerning DGtal dependencies (boost, Qt,...), all the dependencies used to compile your DGtal library must be present to build the DGtalTools.
  
  
