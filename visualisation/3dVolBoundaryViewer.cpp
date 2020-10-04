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
 * @file 3dVolBoundaryViewer.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/11/15
 *
 * A tool file named 3dVolBoundaryViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;


/**
 @page Doc3dVolBoundaryViewer 3dVolBoundaryViewer
 
 @brief  Display the boundary of a volume file by using QGLviewer.
 
 The mode  specifies if you wish to see surface elements (BDRY), the inner
 voxels (INNER) or the outer voxels (OUTER) that touch the boundary.
 
 @b Usage:   3dVolBoundaryViewer -i [input]
 
 @b Allowed @b options @b are :
 
 @code
 
 Positionals:
 1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
 
 Options:
 -h,--help                             Print this help message and exit
 -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
 -m,--thresholdMin INT=0               threshold min (excluded) to define binary shape.
 -M,--thresholdMax INT=255             threshold max (included) to define binary shape.
 --mode TEXT:{INNER,OUTER,BDRY}=INNER  set mode for display: INNER: inner voxels, OUTER: outer voxels, BDRY: surfels
 OUTER: outer voxels, BDRY: surfels
 @endcode
 
 
 @b Example:
 
 
 @code
 3dVolBoundaryViewer  -i $DGtal/examples/samples/lobster.vol -m 60
 @endcode
 
 You should obtain such a result:
 
 @image html res3dVolBoundaryViewer.png "Resulting visualization."
 
 
 @see
 @ref 3dVolBoundaryViewer.cpp
 
 */



int main( int argc, char** argv )
{
  typedef SpaceND<3,int> Space;
  typedef KhalimskySpaceND<3,int> KSpace;
  typedef HyperRectDomain<Space> Domain;
  typedef ImageSelector<Domain, unsigned char>::Type Image;
  typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet;
  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Display the boundary of a volume file by using QGLviewer. The mode specifies if you wish to see surface elements (BDRY), the inner voxels (INNER) or the outer voxels (OUTER) that touch the boundary. \n \t Example: 3dVolBoundaryViewer  -i $DGtal/examples/samples/lobster.vol -m 60");
  std::string inputFileName;
  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  int dicomMin {-1000};
  int dicomMax {3000};
  int thresholdMin {0};
  int thresholdMax {255};
  std::string mode {"INNER"};
  
  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_option("--thresholdMin,-m", thresholdMin, "threshold min (excluded) to define binary shape.", true);
  app.add_option("--thresholdMax,-M", thresholdMax, "threshold max (included) to define binary shape.", true);
#ifdef WITH_ITK
  app.add_option("--dicomMin",dicomMin,"set minimum density threshold on Hounsfield scale", true );
  app.add_option("--dicomMax",dicomMin,"set maximum density threshold on Hounsfield scale", true );
#endif
  app.add_option("--mode", mode,"set mode for display: INNER: inner voxels, OUTER: outer voxels, BDRY: surfels", true )
  -> check(CLI::IsMember({"INNER", "OUTER", "BDRY"}));
  
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  
  
  QApplication application(argc,argv);
  
  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);
  if(extension!="vol" && extension != "p3d" && extension != "pgm3D" && extension != "pgm3d" && extension != "sdp" && extension != "pgm"
#ifdef WITH_ITK
     && extension !="dcm"
#endif
     ){
    trace.info() << "File extension not recognized: "<< extension << std::endl;
    return 0;
  }
  
  if(extension=="vol" || extension=="pgm3d" || extension=="pgm3D"
#ifdef WITH_ITK
     || extension =="dcm"
#endif
     ){
    trace.beginBlock( "Loading image into memory." );
#ifdef WITH_ITK
    int dicomMin = vm["dicomMin"].as<int>();
    int dicomMax = vm["dicomMax"].as<int>();
    typedef DGtal::functors::Rescaling<int ,unsigned char > RescalFCT;
    Image image = extension == "dcm" ? DicomReader< Image,  RescalFCT  >::importDicom( inputFilename,
                                                                                      RescalFCT(dicomMin,
                                                                                                dicomMax,
                                                                                                0, 255) ) :
    GenericReader<Image>::import( inputFilename );
#else
    Image image = GenericReader<Image>::import (inputFileName );
#endif
    trace.info() << "Image loaded: "<<image<< std::endl;
    trace.endBlock();
    
    //! [3dVolBoundaryViewer-KSpace]
    trace.beginBlock( "Construct the Khalimsky space from the image domain." );
    Domain domain = image.domain();
    KSpace ks;
    bool space_ok = ks.init( domain.lowerBound(), domain.upperBound(), true );
    if (!space_ok)
    {
      trace.error() << "Error in the Khamisky space construction."<<std::endl;
      return 2;
    }
    trace.endBlock();
    //! [3dVolBoundaryViewer-KSpace]
    
    //! [3dVolBoundaryViewer-Set3D]
    trace.beginBlock( "Wrapping a digital set around image. " );
    typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
    ThresholdedImage thresholdedImage( image, thresholdMin, thresholdMax );
    trace.endBlock();
    //! [3dVolBoundaryViewer-Set3D]
    
    //! [3dVolBoundaryViewer-ExtractingSurface]
    trace.beginBlock( "Extracting boundary by scanning the space. " );
    typedef KSpace::SurfelSet SurfelSet;
    typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
    typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
    MySurfelAdjacency surfAdj( true ); // interior in all directions.
    MySetOfSurfels theSetOfSurfels( ks, surfAdj );
    Surfaces<KSpace>::sMakeBoundary( theSetOfSurfels.surfelSet(),
                                    ks, thresholdedImage,
                                    domain.lowerBound(),
                                    domain.upperBound() );
    MyDigitalSurface digSurf( theSetOfSurfels );
    trace.info() << "Digital surface has " << digSurf.size() << " surfels."
    << std::endl;
    trace.endBlock();
    //! [3dVolBoundaryViewer-ExtractingSurface]
    
    //! [3dVolBoundaryViewer-ViewingSurface]
    trace.beginBlock( "Displaying everything. " );
    Viewer3D<Space,KSpace> viewer(ks);
    viewer.setWindowTitle("Simple boundary of volume Viewer");
    viewer.show();
    typedef MyDigitalSurface::ConstIterator ConstIterator;
    if ( mode == "BDRY" )
    {
      viewer << SetMode3D(ks.unsigns( *(digSurf.begin()) ).className(), "Basic");
      for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
        viewer << ks.unsigns( *it );
    }else if ( mode == "INNER" )
      for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
        viewer << ks.sCoords( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) );
    else if ( mode == "OUTER" )
      for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
        viewer << ks.sCoords( ks.sIndirectIncident( *it, ks.sOrthDir( *it ) ) );
    else{
      trace.error() << "Warning display mode (" << mode << ") not implemented." << std::endl;
      trace.error() << "The display will be empty." << std::endl;
    }
    viewer << Viewer3D<>::updateDisplay;
    trace.endBlock();
    return application.exec();
  }
  return 0;
}
