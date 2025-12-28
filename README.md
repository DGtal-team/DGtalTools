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
own source code (in complement of DGtal tutorial https://dgtal-team.github.io/doc-nightly/packageTutorials.html).


More Information
----------------
* Related DGtalTools-contrib: https://github.com/DGtal-team/DGtalTools-contrib
* Release 2.1 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17953751.svg)](https://doi.org/10.5281/zenodo.17953751)
* Release 2.0 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17953693.svg)](https://doi.org/10.5281/zenodo.17953693)
* Release 1.4 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11577409.svg)](https://doi.org/10.5281/zenodo.11577409)
* Release 1.3 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7395984.svg)](https://doi.org/10.5281/zenodo.7395984)
* Release 1.2 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4893597.svg)](https://doi.org/10.5281/zenodo.4893597)
* Release 1.1 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4079375.svg)](https://doi.org/10.5281/zenodo.4079375)
* Release 1.0 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2704747.svg)](https://doi.org/10.5281/zenodo.2704747)
* Release 0.9.4.1 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1207828.svg)](https://doi.org/10.5281/zenodo.1207828)
* Release 0.9.4 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1203421.svg)](https://doi.org/10.5281/zenodo.1203421)
* Release 0.9.3 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.290554.svg)](https://doi.org/10.5281/zenodo.290554)
* Release 0.9.2 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.56432.svg)](http://dx.doi.org/10.5281/zenodo.56432)
* Release 0.9.1 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45130.svg)](http://dx.doi.org/10.5281/zenodo.45130)
* Release 0.9 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.31882.svg)](http://dx.doi.org/10.5281/zenodo.31882)
* Release 0.8 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11587.svg)](http://dx.doi.org/10.5281/zenodo.11587)
* Continuous Integration (Linux/MacOS/Windows) ![Build status](https://github.com/DGtal-team/DGtalTools/actions/workflows/buildAndDocumentation.yml/badge.svg)


Documentation
==============

The description and documentation of the tools are available [here](https://dgtal-team.github.io/doctools-nightly). 

Actually the DGTal project is organized as follows:

 - [Converters](https://dgtal-team.github.io/doctools-nightly/converters.html#converters):
   utilities to convert various simple file formats (for instance vol2raw, dicom2vol, mesh2heightfield ...)
   
 - [Estimators](https://dgtal-team.github.io/doctools-nightly/estimators.html#estimators_Doc):
   different geometric estimators (like tangent, curvature 2D/3D...)


 - [Generators](https://dgtal-team.github.io/doctools-nightly/generators.html):
   utilities to generate various contours/shapes


 - [Visualization](https://dgtal-team.github.io/doctools-nightly/visualization.html#visualization_Doc):
   various tools to visualize digital data (set of voxels, vol file, heightmap ... )
   

 - [Volumetric](https://dgtal-team.github.io/doctools-nightly/volumetric.html#volumetric_Doc): 
   tools to manipulate volumetric files (marching cube, sub sampling, thinning)

 - [ImageProcessing](https://dgtal-team.github.io/doctools-nightly/imageProcessing.html): 
   tools to process images (image restoration, image inpainting)


How to build the tools
======================
  - use cmake tool to generate a build script (MakeFile, VS project,..) from the CMakeLists.txt
  - DGtal must be installed in your system. Concerning DGtal dependencies (boost, Qt,...), all the dependencies used to compile your DGtal library must be present to build the DGtalTools.
  
  


Galleries
=========

 - [Converters](https://dgtal-team.github.io/doctools-nightly/converters.html#converters) :
   <center>
   <table>
   <tr>
   <td> <img height=150 src="doc/images/resVol2heightfield.png"> </td>
   <td> <img height=150 src="doc/images/resHeightfield2vol.png"> </td>
   <td> <img height=150 src="doc/images/resImgAddNoise.png"> </td>
   <td> <img height=150 src="doc/images/resVol2slice.png"></td>
   </tr>
   <tr>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__convertertools.html#vol2heightfield_sec">vol2heighfield</a></td>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__convertertools.html#heightfield2vol_sec">heightfield2vol</a></td>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__convertertools.html#imgAddNoise_sec">imgAddNoise</a></td> 
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__convertertools.html#vol2slice_sec">vol2slice</a></td>
   </tr>
   <tr>
   <td> <img height=150 src="doc/images/resMesh2heightfield.png"> 
   <td align=center> <img height=150 src="doc/images/resHeightfield2shading.png"> </td>
   <td> <img height=150 src="doc/images/resVol2sdp.png"></td>
   </tr>
   <tr>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__convertertools.html#mesh2heightfield_sec"> mesh2heightfield</a></td>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__convertertools.html#heightfield2shading_sec"> heightfield2shading</a></td>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__convertertools.html#vol2sdp_sec"> vol2sdp</a></td>
   </tr>
   <tr>
   </table>
   </center>
 
 
 - [Estimators](https://dgtal-team.github.io/doctools-nightly/estimators.html#estimators_Doc) :
   <center>
   <table>
   <tr>
   <td  align="center" colspan="3"><img height=130 src="doc/images/resCurvatureScaleSpaceBCC.png"></td>
   </tr>
   <tr>
   <td colspan="3" align=center >Illustration of <a href="https://dgtal-team.github.io/doctools-nightly/group__estimatortools.html#curvatureScaleSpaceBCC_sec"> curvatureScaleSpaceBCC</a> </td>
   </tr>
   <tr>
   <td align=center ><img height=200 src="doc/images/resGeneric3dNormalEstimatorsVCM.png"> </td>
   <td align=center ><img height=200 src="doc/images/resGeneric3dNormalEstimatorsNoiseVCM.png"> </td>
   <td align=center> <img height=200 src="doc/images/res2dLocalEstimators.png"> </td>
   </tr>
   <tr>
   <td colspan="2"> Illustration of <a href="https://https://dgtal-team.github.io/doctools-nightly/group__estimatortools.html#generic3dNormalEstimators_sec">generic3dNormalEstimators</a> on VCM estimator applied on smooth and noisy shapes.</td>
   <td> <a align=center href="https://dgtal-team.github.io/doctools-nightly/group__estimatortools.html#Doc2dLocalEstimators_sec"> 2dLocalEstimators </a> </td>
   </tr>
   <tr> 
   <td> <img height=200 src="doc/images/res3dCurveTangentEstimator.png"></td>
   <td> <img height=200 src="doc/images/resVol2normalField.png"></td>
   <td> <img height=200 src="doc/images/resLengthEstimators.png"></td>
   </tr>
   <tr> 
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__estimatortools.html#Doc3dCurveTangentEstimator_sec"> 3dCurveTangentEstimator</a> </td>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__estimatortools.html#vol2normalField_sec"> vol2normalField</a> </td>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__estimatortools.html#lengthEstimators_sec"> lengthEstimators</a> </td>
   </tr>
   </table>
   </center>

 - [Generators](https://dgtal-team.github.io/doctools-nightly/generators.html) :
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
   <td colspan="3" align=center> <a href="https://dgtal-team.github.io/doctools-nightly/group__generatorstools.html#shapeGenerator_sec">shapeGenerator</a>  
   </tr>
   </table>
   </center>
   
 - [Visualization](https://dgtal-team.github.io/doctools-nightly/visualization.html#visualization_Doc) :
   <center>
   <table>
   <tr>
   <td><img height=130 src="doc/images/res3dCurvatureViewerNoise.png"></td>
   <td><img height=130 src="doc/images/res3dCurveViewer.png"></td>
   <td><img height=130 src="doc/images/res3dImageViewer.png" ></td>
   <td><img height=130 src="doc/images/res3dVolViewer.png"></td>
   </tr>
   <tr>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3DCurvatureViewer_sec" > 3dCurvatureViewer</a></td>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3dCurveViewer_sec" > 3dCurveViewer</a> </td>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3dImageViewer_sec" > 3dImageViewer</a></td>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3dVolViewer_sec" > 3dVolViewer</a></td>
   </tr>
   <tr>
   <td><img height=130 src="doc/images/resDisplayContours.png"></td>
   <td><img height=130 src="doc/images/resMeshViewer.png"></td>
   <td><img width=300 src="doc/images/res3DSDPViewer.png"> </td>
   <td align="center" colspan="2"><img width=250 src="doc/images/volscope-surface.png"</td>


</tr>
   <tr>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#displayContours"> displayContours</a></td>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#meshViewer"> meshViewer</a></td>
   <td align=center > <a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3DSDPViewer"> 3dSDPViewer</a></td>
   <td  align="center" colspan="2"><a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#DocVolscope_sec">volScope</a>  </td>
   </tr>
   <tr>
   <td><img  width=300 src="doc/images//res3dHeightMapViewer.png" </td>
   <td><img width=300 src="doc/images//res3dCompSurfelData.png" </td>
   <td align="center" colspan="2"><img width=250 src="doc/images//res3dImplicitSurfaceExtractorByThickening.png"</td>
   </tr>

<tr>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3dHeightMapViewer_sec" > 3dHeightMapViewer</a></td>
   <td align=center  ><a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3DCompSurfelData_sec" > 3dCompSurfelData</a></td>
   <td  align="center" colspan="2"><a href="https://dgtal-team.github.io/doctools-nightly/group__visualizationtools.html#Doc3dImplicitSurfaceExtractorByThickening_sec" > 3dImplicitSurfaceExtractorByThickening</a> </td>
   </tr>
<tr>
</tr>
      <tr>

</tr>
   </table>
   </center>





 - [Volumetric](https://dgtal-team.github.io/doctools-nightly/volumetric.html#volumetric_Doc) :
   <center>
   <table>
   <tr>
   <td> <img width=300 src="doc/images/res3dVolMarchingCubes2.png"> </td>
   <td> <img width=300 src="doc/images/resHomotopicThinning3D.png"> </td>
   <td> <img width=300 src="doc/images/resVolSubSample.png"> </td>
   </tr>
   <tr>
   <td align="center"><a href="https://dgtal-team.github.io/doctools-nightly/group__volumetrictools.html#Doc3dVolMarchingCubes_sec"> 3dVolMarchingCubes </a> </td>
   <td align="center"><a href="https://dgtal-team.github.io/doctools-nightly/group__volumetrictools.html#homotopicThinning3D_sec"> homotopicThinning3D </a> </td>
   <td align="center"><a href="https://dgtal-team.github.io/doctools-nightly/group__volumetrictools.html#volSubSample_sec"> volSubSample </a> </td>
   </tr>
   <tr>
   <td align=center > <img width=300 src="doc/images/resVolReSample.png"> </td>
   <td align=center > <img width=300 src="doc/images/resVolTrValues.png"> </td>
   <td align=center > <img width=300 src="doc/images/resVolSubSample.png"> </td>
   </tr>
   <tr>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__volumetrictools.html#volReSample_sec"> volReSample </a> </td>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__volumetrictools.html#volTrValues_sec"> volTrValues </a> </td>
   <td align=center ><a href="https://dgtal-team.github.io/doctools-nightly/group__volumetrictools.html#volSegment_sec"> volSegment </a> </td>
   </tr>
   <tr>
   <td align=center > <img width=300 src="doc/images/resCriticalKernelsThinning3D_select-dmax_skel-1isthmus_persistence-1.png"> </td>
   </tr>
   <tr>
   <td align=center  ><a href="https://dgtal-team.github.io/doctools-nightly/group__volumetrictools.html#criticalKernelsThinning3D_sec"> criticalKernelsThinning3D </a> </td>
   </tr>
   </table>
   </center>

