/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file 3dHeightMapViewer.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/06/2014
 *
 * An example file named 3dHeighMapViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <climits>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/viewers/PolyscopeViewer.h"
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "DGtal/io/Color.h"
#include "DGtal/images/ImageLinearCellEmbedder.h"
#include "DGtal/images/ImageSelector.h"
#include <DGtal/base/ConstAlias.h>
#include <DGtal/kernel/BasicPointFunctors.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/io/colormaps/GradientColorMap.h>
#include <DGtal/io/colormaps/GrayscaleColorMap.h>
#include <DGtal/topology/SetOfSurfels.h>
#include <DGtal/topology/DigitalSurface.h>
#include "DGtal/shapes/TriangulatedSurface.h"
#include "DGtal/shapes/MeshHelpers.h"

using namespace std;
using namespace DGtal;
using namespace Z3i;

#include "CLI11.hpp"



/**
 @page Doc3dHeightMapViewer 3dHeightMapViewer
 
 @brief  Displays 2D image as heightmap by using QGLviewer.
 
 @b Usage:  3dImageViewer [options] input
 
 
 @b Allowed @b options @b are :
 
 @code
 
 Positionals:
 1 TEXT:FILE REQUIRED                  2d input image representing the height map (given as grayscape image cast into 8 bits).
 
 Options:
 -h,--help                             Print this help message and exit
 -i,--input TEXT:FILE REQUIRED         2d input image representing the height map (given as grayscape image cast into 8 bits).
 -s,--scale FLOAT                      set the scale of the maximal level. (default 1.0)
 -c,--colorMap                         define the heightmap color with a pre-defined colormap (GradientColorMap)
 -t,--colorTextureImage TEXT           define the heightmap color from a given color image (32 bits image).
 
 @endcode
 
 
 @b Example:
 
 @code
 $ visualisation/3dHeightMapViewer ${DGtal}/examples/samples/church.pgm -s 0.2
 @endcode
 
 You should obtain such a result:
 
 @image html res3dHeightMapViewer.png "resulting heightmap visualisation."
 
 
 @see
 @ref 3dHeightMapViewer.cpp
 
 */

// Defining a Helper to get the 3D point functor from an 2DImage
template<typename TImage2D, typename TPoint3D >
struct Image3DPredicatFrom2DImage{
    typedef  TPoint3D Point3D;
    /**
     *  Construct the predicat given a 2D Image
     **/
    Image3DPredicatFrom2DImage(DGtal::ConstAlias<TImage2D> anImage, double aScale):myImageRef(anImage),
    myScale(aScale){
    }
    inline
    bool operator()(const Point3D &aPoint)  const {
        functors::Projector<SpaceND<2, typename TImage2D::Integer> > projXY;
        return  (*myImageRef)(projXY(aPoint))*myScale >= aPoint[2];
    }
    CountedConstPtrOrConstPtr<TImage2D> myImageRef;
    double myScale;
};




typedef ImageSelector < Domain, int>::Type Image;



int main( int argc, char** argv )
{
    
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    std::string inputFileName;
    double scale {1.0};
    bool colorMap {false};
    std::string colorTextureImage;
    app.description("Displays 2D image as heightmap by using QGLviewer.\n Exemple of use:  visualisation/3dHeightMapViewer  ${DGtal}/examples/samples/church.pgm -s 0.2");
    
    app.add_option("-i,--input,1", inputFileName, "2d input image representing the height map (given as grayscape image cast into 8 bits)." )
    ->required()
    ->check(CLI::ExistingFile);
    app.add_option("--scale,-s",scale,  "set the scale of the maximal level. (default 1.0)");
    app.add_flag("--colorMap,-c", colorMap, "define the heightmap color with a pre-defined colormap (GradientColorMap)");
    app.add_option("--colorTextureImage,-t", colorTextureImage,  "define the heightmap color from a given color image (32 bits image).");
    
    
    
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------
    
    
    typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2DG ;
    typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned int> Image2DCol ;
    
    Image2DG image = GenericReader<Image2DG>::import( inputFileName );
    Image2DCol imageTexture(image.domain());
    int maxHeight = (int)(std::numeric_limits<Image2DG::Value>::max()*scale);
    trace.info()<< "Max height from scale:" << maxHeight << std::endl;
    GradientColorMap<Image2DG::Value,CMAP_JET>  gradientShade( 0, std::numeric_limits<Image2DG::Value>::max());
    GrayscaleColorMap<Image2DG::Value>  grayShade(0, std::numeric_limits<Image2DG::Value>::max());

    if(colorTextureImage != ""){
        imageTexture =  GenericReader<Image2DCol>::import( colorTextureImage );
    }
    
    Z2i::RealPoint plow (image.domain().lowerBound()[0]-0.5,
                         image.domain().lowerBound()[1]-0.5);
    
    Z2i::RealPoint pup (image.domain().upperBound()[0]+0.5,
                        image.domain().upperBound()[1]+0.5);
    
    PolyscopeViewer viewer;
    bool intAdjacency = true;
    
    typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
    typedef KSpace::SurfelSet SurfelSet;
    typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
    typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
    
    typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
    MySurfelAdjacency surfAdj( intAdjacency ); // interior in all directions.
    
    
    KSpace K;
    K.init(Z3i::Point(0,0,0),Z3i::Point(image.domain().upperBound()[0], image.domain().upperBound()[1], maxHeight+1), true);
    SurfelSet boundVect;
    Image3DPredicatFrom2DImage<Image2DG, Z3i::Point> image3Dpredicate(image, scale);
    trace.info() << "Constructing boundary... ";
    MySetOfSurfels theSetOfSurfels( K, surfAdj );
    
    Surfaces<KSpace>::sMakeBoundary( theSetOfSurfels.surfelSet(),
                                    K, image3Dpredicate, Z3i::Point(0,0,0),
                                    Z3i::Point(image.domain().upperBound()[0], image.domain().upperBound()[1], maxHeight+1) );
    trace.info() << "[done]"<< std::endl;
    
    MyDigitalSurface digSurf( theSetOfSurfels );
    trace.info() << "Digital surface has " << digSurf.size() << " surfels."
    << std::endl;
    
    trace.beginBlock( "Making triangulated surface. " );
    typedef CanonicEmbedder< Space >                                  TrivialEmbedder;
    typedef ImageLinearCellEmbedder< KSpace,Image , TrivialEmbedder > CellEmbedder;
    typedef CellEmbedder::Value                                       RealPoint;
    typedef TriangulatedSurface< RealPoint >                          TriMesh;
    typedef Mesh< RealPoint >                                         ViewMesh;
    typedef std::map< MyDigitalSurface::Vertex, TriMesh::Index >      VertexMap;
    TriMesh         trimesh;
    ViewMesh        viewmesh;
    TrivialEmbedder trivialEmbedder;
    CellEmbedder    cellEmbedder;
    Image::Domain d (Z3i::Point(0,0,0),Z3i::Point(image.domain().upperBound()[0], image.domain().upperBound()[1], maxHeight+1));
    Image imageE(d);
    for (auto p : imageE.domain()){
        imageE.setValue(p, image3Dpredicate(p));
    }
    cellEmbedder.init(K, imageE, trivialEmbedder, 0);
    VertexMap vmap;
    MeshHelpers::digitalSurface2DualTriangulatedSurface
    ( digSurf, cellEmbedder, trimesh, vmap );
    MeshHelpers::triangulatedSurface2Mesh( trimesh, viewmesh );
    trace.info() << "Mesh has " << viewmesh.nbVertex()
    << " vertices and " << viewmesh.nbFaces() << " faces." << std::endl;
    trace.endBlock();
    for (unsigned int i = 0; i < viewmesh.nbFaces(); i++){
        auto b =  viewmesh.getFaceBarycenter(i);
        auto bp = Z2i::Point(b[0], b[1]);
        auto val = image(bp);
        if (image.domain().isInside(bp)){
            if(colorMap){
                viewmesh.setFaceColor(i,gradientShade(val));
            }else if (colorTextureImage != "") {
                viewmesh.setFaceColor(i, DGtal::Color(imageTexture(bp)));
            }else{
                viewmesh.setFaceColor(i, grayShade(val));
            }
                        
        }
    }
    viewer.drawColor(Color(150,0,0,254));
    viewer << viewmesh;
    
    viewer.show();
    return 0;
}

