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


More Information
----------------
* Related DGtalTools-contrib: https://github.com/DGtal-team/DGtalTools-contrib
* Release 0.9.1 DOI [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45130.svg)](http://dx.doi.org/10.5281/zenodo.45130)
* Release 0.9 DOI [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.31882.svg)](http://dx.doi.org/10.5281/zenodo.31882)
* Release 0.8 DOI [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11587.svg)](http://dx.doi.org/10.5281/zenodo.11587)
* [![Build Status](https://api.travis-ci.org/DGtal-team/DGtalTools.svg?branch=master)](https://travis-ci.org/DGtal-team/DGtalTools)


Documentation
==============

The description and documentation of the tools are available [here] (http://dgtal.org/doc/tools/nightly/). 

Actually the DGTal project is organized as follows:

 - [Converters](http://dgtal.org/doc/tools/nightly/converters.html#converters)
    Utilities to convert various simple file formats (for instance vol2raw, dicom2vol, mesh2heightfield ...). 
    <center>
    <table>
    <tr>
    <td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/6706730/9bac9720-cd60-11e4-9819-81e536b21e97.gif"> </td>
    <td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/15342702/d9999fd8-1c97-11e6-8f6a-6bbd114a5641.png"> </td>
    <td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/6914017/ccf433e4-d786-11e4-997b-f513f07f56f3.gif"> </td>
    <td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/6914104/6d25be78-d787-11e4-9433-2a834fc4a0af.png"> </td>
    </tr>
    <tr>
    <td> <a href="http://dgtal.org/doc/tools/nightly/vol2heightfield.html">vol2heightmap</a></td>
    <td>[heightfield2vol](http://dgtal.org/doc/tools/nightly/heightfield2vol.html)</td>
    <td>[imgAddNoise](http://dgtal.org/doc/tools/nightly/imgAddNoise.html)</td> 
    <td>[volAddNoise](http://dgtal.org/doc/tools/nightly/volAddNoise.html)</td>
    </tr>
    <tr>
    <td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/7106643/c09cb8d4-e148-11e4-8653-2d5bac5dc3c5.png"> 
    <td align=center> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/10265156/64aaad64-6a24-11e5-9773-628c661fa76c.png"> </td>
    <td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/15342501/d80656d6-1c95-11e6-8461-2ec25ba5a864.png"></td>
    <td> <img height=150 src="https://cloud.githubusercontent.com/assets/772865/15342675/b17e505c-1c97-11e6-85c3-663550b3b39a.png"></td>
    </tr>
    <tr>
    <td>[mesh2heightfield](http://dgtal.org/doc/tools/nightly/mesh2heightfield.html)</td>
    <td>[heightfield2shading](http://dgtal.org/doc/tools/nightly/heightfield2shading.html)</td>
    <td>[vol2sdp](http://dgtal.org/doc/tools/nightly/vol2sdp.html)</td>
    <td>[vol2slice](http://dgtal.org/doc/tools/nightly/vol2slice.html)</td>
    </tr>
    <tr>
    </table>
    </center>
    
 - [Estimators](http://dgtal.org/doc/tools/nightly/estimators.html#estimators_Doc)
    Different geometric estimators (like tangent, curvature 2D/3D...)
    <table>
    <tr>
    <td colspan="3"><img height=130 src="https://cloud.githubusercontent.com/assets/772865/2646108/f515b0a2-bf39-11e3-96f8-c7606173f43b.png"></td>
    </tr>
    <tr>
    <td colspan="3">Illustration of curvatureScaleSpaceBCC </td>
    </tr>
    <tr>
    <td align=center ><img height=200 src="https://cloud.githubusercontent.com/assets/793707/2996392/d3ee9e58-dced-11e3-98a0-72233927aaf6.jpg"> </td>
    <td align=center ><img height=200 src="https://cloud.githubusercontent.com/assets/772865/15353086/8637153c-1ce7-11e6-9b64-45baaa954d98.png"> </td>
    <td align=center> <img height=200 src="https://cloud.githubusercontent.com/assets/772865/15352637/6cfde4a8-1ce5-11e6-8d36-a2a3be4744d1.png"> </td>
    </tr>
    <tr>
    <td colspan="2"> Illustration of generic3dNormalEstimators on VCM estimator applied on smooth and noisy shapes.</td>
    <td> 2dLocalEstimators </td>
    </tr>
    <tr> 
    <td> <img height=200 src="https://cloud.githubusercontent.com/assets/772865/15351379/5d2cd562-1cdf-11e6-81b6-8171b0e7f87f.png"></td>
    <td> <img height=200 src="https://cloud.githubusercontent.com/assets/772865/15351501/dbb860d6-1cdf-11e6-8d9f-df3356bb969c.png"></td>
    <td> <img height=200 src="https://cloud.githubusercontent.com/assets/772865/15352775/29e2fad6-1ce6-11e6-958c-34ca77175a90.png"></td>
    </tr>
    <tr> 
    <td> [3dCurveTangentEstimator](http://dgtal.org/doc/tools/nightly/Doc3dCurveTangentEstimator.html) </td>
    <td> [vol2normalField](http://dgtal.org/doc/tools/nightly/vol2normalField.html) </td>
    <td> [lengthEstimators](http://dgtal.org/doc/tools/nightly/lengthEstimators.html) </td>
    </tr>
    </table>


 - [Generators](http://dgtal.org/doc/tools/nightly/generators.html)
    Utilities to generate various contours/shapes:
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
    <td colspan=3 > Illustration of the Generators tools </td>
    </tr>
    </table>
    </center>


 - [Visualization](http://dgtal.org/doc/tools/nightly/visualization.html#visualization_Doc)
    various tools to visualize digital data (set of voxels, vol file, heightmap ... )
    <center>
    <table>
    <tr>
    <td><img height=130 src="https://f.cloud.github.com/assets/772865/684607/450c064e-da00-11e2-8830-76eb90a5efd7.png"></td>
    <td><img height=130 src="https://f.cloud.github.com/assets/772865/685853/d96a5252-da44-11e2-9872-7f0160be8f5d.png"></td>
    <td><img height=130 src="https://f.cloud.github.com/assets/772865/684569/59a2f6fa-d9fe-11e2-84ba-a48842f4aafb.png" ></td>
    <td><img height=130 src="https://f.cloud.github.com/assets/772865/684590/778bea9a-d9ff-11e2-8e04-6e3e8a39ae3c.png"></td>
    </tr>
    <tr>
    <td>[3dCurvatureViewer](http://dgtal.org/doc/tools/nightly/Doc3DCurvatureViewer.html)</td>
    <td>[3dCurveViewer](http://dgtal.org/doc/tools/nightly/Doc3dCurveViewer.html) </td>
    <td>[3dImageViewer](http://dgtal.org/doc/tools/nightly/Doc3dImageViewer.html)</td>
    <td>[3dVolViewer](http://dgtal.org/doc/tools/nightly/Doc3dVolViewer.html)</td>
    </tr>

    <tr>
    <td><img height=130 src="https://f.cloud.github.com/assets/772865/684598/c3adcf4c-d9ff-11e2-8c3f-e67c8abd0c76.png"></td>
    <td><img height=130 src="https://f.cloud.github.com/assets/772865/684622/d698405a-da00-11e2-8aa0-19212a58ce23.png"></td>
    <td><img width=300 src="https://cloud.githubusercontent.com/assets/772865/2720141/6c42a0e0-c56b-11e3-8328-a6d88242f21e.png"> </td>
    <td><img  width=300 src="https://cloud.githubusercontent.com/assets/772865/10269635/8678ca02-6add-11e5-9d83-fbcf3608612f.png"></td>
    </tr>
    <tr>
    <td>[displayContours](http://dgtal.org/doc/tools/nightly/displayContours.html)</td>
    <td>[meshViewer](http://dgtal.org/doc/tools/nightly/meshViewer.html)</td>
    <td>[3dSDPViewer](http://dgtal.org/doc/tools/nightly/Doc3DSDPViewer.html)</td>
    <td>[sliceViewer](http://dgtal.org/doc/tools/nightly/sliceViewer.html)</td>
    </tr>

    <tr>
    <td><img  width=300 src="https://cloud.githubusercontent.com/assets/772865/3486505/6edb2144-043e-11e4-81c4-2c20f272a119.png" </td>
    <td><img width=300 src="https://cloud.githubusercontent.com/assets/772865/4118303/5f53f61c-329f-11e4-9629-23c53afd9eff.png" </td>
    <td align="center" colspan="2"><img width=250 src="https://cloud.githubusercontent.com/assets/772865/12538898/00459c06-c2e6-11e5-91f6-d14494ab05da.png"</td>
    </tr>
    <tr>
    <td>[3dHeightMapViewer](http://dgtal.org/doc/tools/nightly/Doc3dHeightMapViewer.html)</td>
    <td>[3dCompSurfelData](http://dgtal.org/doc/tools/nightly/CompSurfelData.html)</td>
    <td  align="center" colspan="2">[3dImplicitSurfaceExtractorByThickening](http://dgtal.org/doc/tools/nightly/Doc3dImplicitSurfaceExtractorByThickening.html) </td>
    </tr>
    </table>
    </center>

   


 - [Volumetric](http://dgtal.org/doc/tools/nightly/volumetric.html#volumetric_Doc)
   Tools to manipulate volumetric files (marching cube, sub sampling, thinning).

    <center>
    <table>
    <tr>
    <td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/15342874/d79fc22e-1c98-11e6-922b-48f27f5cf586.png"> </td>
    <td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/15342814/828267ec-1c98-11e6-9080-a210bdd42b35.png"> </td>
    <td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/15342834/a6736bb0-1c98-11e6-9a42-9dc5a198800f.png"> </td>
    </tr>
    <tr>
    <td>[3dVolMarchingCubes](http://dgtal.org/doc/tools/nightly/Doc3dVolMarchingCubes.html)</td>
    <td>[homotopicThinning3D](http://dgtal.org/doc/tools/nightly/homotopicThinning3D.html)</td>
    <td>[volSubSample](http://dgtal.org/doc/tools/nightly/volSubSample.html)</td>
    </tr>
    <tr>
    <td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/15342933/2da32a80-1c99-11e6-9600-16a68f145bdd.png"> </td>
    <td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/15342979/8e2d2d4c-1c99-11e6-88b4-5835adb04eff.png"> </td>
    <td> <img width=300 src="https://cloud.githubusercontent.com/assets/772865/15343002/b2fdf4c6-1c99-11e6-986b-9f56299b4a16.png"> </td>
    </tr>
    <tr>
    <td>[volReSample](http://dgtal.org/doc/tools/nightly/volReSample.html)</td>
    <td>[volTrValues](http://dgtal.org/doc/tools/nightly/volTrValues.html)</td>
    <td>[volSegment](http://dgtal.org/doc/tools/nightly/volSegment.html)</td>
    </tr>
    </table>
    </center>



How to build the tools
======================
  - use cmake tool to generate a build script (MakeFile, VS project,..) from the CMakeLists.txt
  - DGtal must be installed in your system. Concerning DGtal dependencies (boost, Qt,...), all the dependencies used to compile your DGtal library must be present to build the DGtalTools.
  
  
