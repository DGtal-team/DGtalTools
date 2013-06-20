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

 - [Estimators](#estimators)

 - [ShapeGenerator](#shapeGenerator)

 - [Visualization](#visualizatoin)

 - [Volumetric](#volumetric)


Converters
----------

Utilities to convert various simple file formats:

  - freeman2sdp: convert freeman chain towards a Sequence of Discrete Points.
  -  pgm2freeman: to extract a freeman chain contour from a grayscale image.
  - raw2vol and vol2raw: transform 3D volumes files from (resp. to) raw to vol.
  - ofs2off: convert OFS mesh format towards a OFF variant. 


Estimators
----------
 
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
  - estimatorComparator: program to perform comparison of local quantity estimators (to be updated for current DGtal version).
  - vol2normalField: compute the normal vector field of a given vol file .


ShapeGenerator    
--------------
  - shapeGenerator: generate multigrid shape
  - contourGenerator: generate multigrid shape contours

Visualization
-------------
  - 3dCurvatureViewer: permits to compute and visualize mean or gaussian curvature of binary shapes.
  - 3dImageViewer: tools to display 3d slice images with QGLViewer.
  - 3dVolViewer: volume file (.vol and .pgm3d) viewer with QGLViewer.
  - displayContours: display discrete contours from various format (.fc (freemanchain), .sdp).
  - meshViewer: display 3D mesh from OFS or OFF format. 
  - patternTriangulation: a new tool that draws with Board2D the convex hull, the closest-point Delaunay triangulation or the farthest-point Delaunay triangulation of a pattern.
 
  
Here are some illustration of such tools:
<center>
<table>
<tr>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684607/450c064e-da00-11e2-8830-76eb90a5efd7.png"></td>
<td><img width=130 src="https://f.cloud.github.com/assets/772865/684569/59a2f6fa-d9fe-11e2-84ba-a48842f4aafb.png" ></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684590/778bea9a-d9ff-11e2-8e04-6e3e8a39ae3c.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684598/c3adcf4c-d9ff-11e2-8c3f-e67c8abd0c76.png"></td>
<td><img height=130 src="https://f.cloud.github.com/assets/772865/684622/d698405a-da00-11e2-8aa0-19212a58ce23.png"></td>
</tr>
<td>3dCurvatureViewer</td>
<td>3dImageViewer</td>
<td>3dVolViewer</td>
<td>displayContours</td>
<td>meshViewer</td>
</table>
<center>
Volumetric
----------
  - 3dVolMarchingCubes: marching cubes form a Vol file	
  - homotopicThinning3D: ultimate skeleton from vol file
  - volAddBorder: add a 1 voxel boundary with value 0 to a vol file.
  - volCComponentCounter: a simple program to count the number of connected components in a 3D image.
  - volSubSample: sub sample a vol file (division by 2 in each direction).
  


How to build the tools
======================
  - use cmake tool to generate a build script (MakeFile, VS project,..) from the CMakeLists.txt
  - DGtal must be installed in your system. Concerning DGtal dependencies (boost, Qt,...), all the dependencies used to compile your DGtal library must be present to build the DGtalTools.
   
