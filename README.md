
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
own source code (in complement of DGtal tutorial
http://libdgtal.org/doc/stable/tutorials.html).


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
  - 3dVolViewer: volume file (.vol and .pgm3d) viewer with QGLViewer.
  - 3dImageViewer: tools to display 3d slice images with QGLViewer.
  - displayContours: display discrete contours from various format (.fc (freemanchain), .sdp).
  - meshViewer: display 3D mesh from OFS or OFF format. 
  - 3dCurvatureViewer: permits to compute and visualize mean or gaussian curvature of binary shapes.


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
   
