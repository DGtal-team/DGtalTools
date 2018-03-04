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
 * @file volBoundary2obj.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 *
 * @date 2013/11/15
 *
 * A tool file named volBoundary2obj.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <set>
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
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"


using namespace std;
using namespace DGtal;
//using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

/**
   @page volBoundary2obj volBoundary2obj
   @brief  Extracts digital points from 3d vol files.

   @b Usage: volBoundary2obj [input] [output]

   @b Allowed @b options @b are:

   @code
   -h [ --help ]                    display this message
   -i [ --input ] arg               vol file (.vol, .longvol .p3d, .pgm3d and if
                                    WITH_ITK is selected: dicom, dcm, mha, mhd).
                                    For longvol, dicom, dcm, mha or mhd formats,
                                    the input values are linearly scaled between
                                    0 and 255.
   -o [ --output ] arg              output obj file (.obj)
   -m [ --thresholdMin ] arg (=0)   threshold min (excluded) to define binary 
   shape
   -M [ --thresholdMax ] arg (=255) threshold max (included) to define binary 
   shape
   --rescaleInputMin arg (=0)       min value used to rescale the input 
                                    intensity (to avoid basic cast into 8  bits 
                                    image).
   --rescaleInputMax arg (=255)     max value used to rescale the input 
                                    intensity (to avoid basic cast into 8 bits 
                                    image).

   --mode arg (=BDRY)               set mode for display: INNER: inner voxels,
   OUTER: outer voxels, BDRY: surfels, CLOSURE:
   surfels with linels and pointels.
   -n [ --normalization ]           Normalization so that the geometry fits in 
   [-1/2,1/2]^3
   @endcode

   @b Example:
   @code 
   $ volBoundary2obj -i $DGtal/examples/samples/lobster.vol -m 80 -o out.obj
   @endcode

   You should obtain such a visualization:
   @image html resVolBoundary2obj.png "resulting visualisation."

   @see
   @ref volBoundary2obj.cpp

*/


int main( int argc, char** argv )
{
  typedef SpaceND<3,int> Space;
  typedef KhalimskySpaceND<3,int> KSpace;
  typedef HyperRectDomain<Space> Domain;
  typedef ImageSelector<Domain, unsigned char>::Type Image;
  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ("output,o", po::value<std::string>(), "output obj file (.obj)" )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
    ("rescaleInputMin", po::value<DGtal::int64_t>()->default_value(0), "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).")
    ("rescaleInputMax", po::value<DGtal::int64_t>()->default_value(255), "max value used to rescale the input intensity (to avoid basic cast into 8 bits image).")
    ("mode",  po::value<std::string>()->default_value("BDRY"), "set mode for display: INNER: inner voxels, OUTER: outer voxels, BDRY: surfels (default), CLOSURE: surfels with linels and pointels.")
    ("normalization,n", "Normalization so that the geometry fits in [-1/2,1/2]^3") ;

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
      std::cout << "Usage: " << argv[0] << " -i [input] -o [output]\n"
                << "Export the boundary of a volume file to OBJ format. The mode specifies if you wish to see surface elements (BDRY), the inner voxels (INNER) or the outer voxels (OUTER) that touch the boundary."<< endl
                << general_opt << "\n";
      return 0;
    }

  if(! vm.count("input"))
    {
      trace.error() << " The file name was defined" << endl;
      return 0;
    }

  if(! vm.count("output"))
    {
      trace.error() << " The output filename was defined" << endl;
      return 0;
    }

  string inputFilename = vm["input"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();
  string mode = vm["mode"].as<string>();
  bool normalization = false;
  if  (vm.count("normalization"))
    normalization = true;
  
  DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
  DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();

  trace.beginBlock( "Loading file.." );
  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image image =  GenericReader< Image >::importWithValueFunctor( inputFilename,RescalFCT(rescaleInputMin,
                                                                                         rescaleInputMax,
                                                                                         0, 255) );  
  trace.info() << "Image loaded: "<<image<< std::endl;
  trace.endBlock();

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

  trace.beginBlock( "Wrapping a digital set around image. " );
  typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
  ThresholdedImage thresholdedImage( image, thresholdMin, thresholdMax );
  trace.endBlock();

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

  trace.beginBlock( "Exporting everything." );
  Board3D<Space,KSpace> board(ks);

  board << SetMode3D(  ks.unsigns( *digSurf.begin() ).className(), "Basic" );
  board << SetMode3D(  (*digSurf.begin()).className(), "Basic" );

  typedef MyDigitalSurface::ConstIterator ConstIterator;
  if ( mode == "BDRY" )
    for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
    { 
      //board << ks.unsigns( *it );
      bool indirectFilled = thresholdedImage(ks.sCoords( ks.sIndirectIncident( *it, ks.sOrthDir( *it ) ) ));
      bool directFilled = thresholdedImage(ks.sCoords( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) ));
      Z3i::RealPoint sc = board.embedKS(*it);
      Z3i::RealPoint n (0,0,0);
      auto ip = ks.sCoords( ks.sIndirectIncident( *it, ks.sOrthDir( *it ) ) );  
      if (!indirectFilled)
      {
        ip = ks.sCoords( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) );
      }
      n = sc - Z3i::RealPoint(static_cast<Z3i::RealPoint::Component>(ip[0]),
                              static_cast<Z3i::RealPoint::Component>(ip[1]),
                              static_cast<Z3i::RealPoint::Component>(ip[2]));
      bool xodd = ks.sIsOpen( *it, 0 );
      bool yodd = ks.sIsOpen( *it, 1 );
      bool zodd = ks.sIsOpen( *it, 2 );      
      board.addQuadFromSurfelCenterWithNormal
        ( Z3i::RealPoint( sc[0]+(xodd? 0:0.5 ), sc[1]+(yodd? 0:0.5 ), sc[2]+(zodd? 0:0.5 ) ),
          ! xodd, ! yodd, ! zodd, n, true,  true, false );
    }
  else if ( mode == "INNER" )
    for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
      board << ks.sCoords( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) );
  else if ( mode == "OUTER" )
    for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
      board << ks.sCoords( ks.sIndirectIncident( *it, ks.sOrthDir( *it ) ) );
  else  if (mode == "CLOSURE")
  {
    std::set<KSpace::Cell> container;
      for ( ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it )
        {
          container.insert( ks.unsigns( *it ) );
          KSpace::SCells oneNeig = ks.sLowerIncident(*it);
          //Processing linels
          for(KSpace::SCells::ConstIterator itt = oneNeig.begin(), ittend = oneNeig.end(); itt != ittend; ++itt)
            {
              container.insert( ks.unsigns( *itt) );
              KSpace::SCells oneNeig2 = ks.sLowerIncident(*itt);
              //Processing pointels
              for(KSpace::SCells::ConstIterator ittt = oneNeig2.begin(), itttend = oneNeig2.end(); ittt != itttend; ++ittt)
                container.insert( ks.unsigns(*ittt) );
            }
        }
      trace.info()<< "Exporting "<< container.size() << " cells"<<std::endl;
      for(auto cell: container)
        board << cell;
    }
    
  string outputFilename = vm["output"].as<std::string>();

  board.saveOBJ(outputFilename, normalization);
  trace.endBlock();

  return 0;
}
