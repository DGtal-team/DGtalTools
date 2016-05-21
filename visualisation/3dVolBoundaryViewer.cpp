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
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

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
  -h [ --help ]                    display this message
  -i [ --input ] arg               vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm
                                   (with 3 dims)) file or sdp (sequence of 
                                   discrete points)
  -m [ --thresholdMin ] arg (=0)   threshold min (excluded) to define binary 
                                   shape
  -M [ --thresholdMax ] arg (=255) threshold max (included) to define binary 
                                   shape
  --dicomMin arg (=-1000)          set minimum density threshold on Hounsfield 
                                   scale
  --dicomMax arg (=3000)           set maximum density threshold on Hounsfield 
                                   scale
  --mode arg (=INNER)              set mode for display: INNER: inner voxels, 
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

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef SpaceND<3,int> Space;
  typedef KhalimskySpaceND<3,int> KSpace;
  typedef HyperRectDomain<Space> Domain;
  typedef ImageSelector<Domain, unsigned char>::Type Image;
  typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet;
  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
#ifdef WITH_ITK
    ("dicomMin", po::value<int>()->default_value(-1000), "set minimum density threshold on Hounsfield scale")
    ("dicomMax", po::value<int>()->default_value(3000), "set maximum density threshold on Hounsfield scale")
#endif
    ("mode",  po::value<std::string>()->default_value("INNER"), "set mode for display: INNER: inner voxels, OUTER: outer voxels, BDRY: surfels") ;

  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " -i [input]\n"
                << "Display the boundary of a volume file by using QGLviewer. The mode specifies if you wish to see surface elements (BDRY), the inner voxels (INNER) or the outer voxels (OUTER) that touch the boundary."<< endl
                << general_opt << "\n";
      return 0;
    }

  if(! vm.count("input"))
    {
      trace.error() << " The file name was defined" << endl;
      return 0;
    }
  string inputFilename = vm["input"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();
  string mode = vm["mode"].as<string>();

  QApplication application(argc,argv);

  string extension = inputFilename.substr(inputFilename.find_last_of(".") + 1);
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
    Image image = GenericReader<Image>::import (inputFilename );
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
    if ( mode == "BDRY" ){
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
