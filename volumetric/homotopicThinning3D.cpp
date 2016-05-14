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
 * @file homotopicThinning3D.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/04
 *
 * Apply an the homotopic thinning from a volumetric volume file.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <QImageReader>

#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/Color.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/PointListReader.h"

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/images/ImageSelector.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////

/**
 @page homotopicThinning3D homotopicThinning3D
 
 @brief Applies an homotopic thinning of a 3d image file (vol,longvol,pgm3d...) with 3D viewer.

 @b Usage: homotopicThinning3d [options] --input <3dImageFileName>  {vol,longvol,pgm3d...}


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]              display this message.
  -i [ --input ] arg         Input volumetric file (.vol, .pgm3d or p3d)
  -m [ --min ] arg (=0)      Minimum (excluded) value for threshold.
  -M [ --max ] arg (=255)    Maximum (included) value for threshold.
  -e [ --exportSDP ] arg     Export the resulting set of points in a simple 
                             (sequence of discrete point (sdp)).
  --fixedPoints arg          defines the coordinates of points which should not
                             be removed.
  -s [ --fixedPointSDP ] arg use fixed points from a file.

 @endcode

 @b Example: 
 
 Usage by forcing point to be left by the thinning: 
 
@code
 $  homotopicThinning3D --input ${DGtal}/examples/samples/Al.100.vol  --fixedPoints 56 35 5  56 61 5  57 91 38  58 8 38  45 50 97 
 @endcode


 You should obtain such a result:
 @image html resHomotopicThinning3D.png "Resulting visualization."
 
 @see
 @ref homotopicThinning3D.cpp

 */

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( std::string param )
{
    trace.error() <<" Parameter: "<<param<<" is required..";
    trace.info() <<std::endl;
    exit ( 1 );
}



int main( int argc, char** argv )
{

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input volumetric file (.vol, .pgm3d or p3d)" )
    ( "min,m", po::value<int>()->default_value( 0 ), "Minimum (excluded) value for threshold.")
    ( "max,M", po::value<int>()->default_value( 255 ), "Maximum (included) value for threshold.")
    ( "exportSDP,e", po::value<std::string>(), "Export the resulting set of points in a simple (sequence of discrete point (sdp)).")
    ("fixedPoints", po::value<std::vector <int> >()->multitoken(), "defines the coordinates of points which should not be removed." )
    ( "fixedPointSDP,s", po::value<std::string>(), "use fixed points from a file.");


  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if ( !parseOK || vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Illustration of homotopic thinning of a 3d image file (vol,longvol,pgm3d...) with 3D viewer."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\thomotopicThinning3d [options] --input <3dImageFileName>  {vol,longvol,pgm3d...} "<<std::endl
                   << general_opt << "\n"
                   << " Usage by forcing point to be left by the thinning: \n"
                   << "homotopicThinning3D --input ${DGtal}/examples/samples/Al.100.vol  --fixedPoints 56 35 5  56 61 5  57 91 38  58 8 38  45 50 97 \n";




      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();

  
  typedef ImageSelector < Z3i::Domain, unsigned char>::Type Image;
  Image image = GenericReader<Image>::import ( filename );

  trace.beginBlock("DT Computation");
  typedef functors::IntervalForegroundPredicate<Image> Predicate;
  Predicate aPredicate(image, vm[ "min" ].as<int>(), vm[ "max" ].as<int>() );

  const Z3i::L2Metric aMetric{};
  DistanceTransformation<Z3i::Space, Predicate , Z3i::L2Metric> dt(image.domain(), aPredicate, aMetric );
  trace.endBlock();
  trace.info() <<image<<std::endl;

  // Domain creation from two bounding points.

  trace.beginBlock("Constructing Set");
  DigitalSet shape_set( image.domain() );
  DigitalSet fixedSet( image.domain() );

  // Get the optional fixed points
  if( vm.count("fixedPoints")){
    std::vector<int> vectC = vm["fixedPoints"].as<std::vector<int> >();
    if(vectC.size()%3==0){
      for( unsigned int i=0; i < vectC.size()-2; i=i+3)
        {
          Z3i::Point pt(vectC.at(i), vectC.at(i+1), vectC.at(i+2));
          fixedSet.insertNew(pt);
        }
    }else{
      trace.error()<< " The coordinates should be 3d coordinates, ignoring fixedPoints option." << std::endl;
    }
  }
  if(vm.count("fixedPointSDP"))
    {
      std::vector<Z3i::Point> vPt = PointListReader<Z3i::Point>::getPointsFromFile(vm["fixedPointSDP"].as<std::string>());
      for( auto &p: vPt)
        {
          fixedSet.insert(p);
        }
    }
  
  SetFromImage<DigitalSet>::append<Image>(shape_set, image,
                                          vm[ "min" ].as<int>(), vm[ "max" ].as<int>() );
  trace.info() << shape_set<<std::endl;
  trace.endBlock();




  trace.beginBlock("Computing skeleton");
  // (6,18), (18,6), (26,6) seem ok.
  // (6,26) gives sometimes weird results (but perhaps ok !).
  Object26_6 shape( dt26_6, shape_set );
  int nb_simple=0;
  int layer = 1;
  std::queue<DigitalSet::Iterator> Q;
  do
    {
      trace.info() << "Layer: "<< layer << std::endl;
      int nb=0;
      DigitalSet & S = shape.pointSet();

      trace.progressBar(0, (double)S.size());
      for ( DigitalSet::Iterator it = S.begin(); it != S.end(); ++it )
        {
          if ( nb % 100 == 0 ) trace.progressBar((double)nb, (double)S.size());
          nb++;
          if (dt( *it ) <= layer)
            {
              if ( shape.isSimple( *it ) )
                Q.push( it );
            }
        }
      trace.progressBar( (double)S.size(), (double)S.size() );
      nb_simple = 0;
      while ( ! Q.empty() )
        {
          DigitalSet::Iterator it = Q.front();
          Q.pop();
          if ( shape.isSimple( *it ) && fixedSet.find(*it) == fixedSet.end() )
            {
              S.erase( *it );
              ++nb_simple;
            }
        }
      trace.info() << "Nb simple points : "<<nb_simple<< " " << std::endl;
      ++layer;
    }
  while ( nb_simple != 0 );
  trace.endBlock();

  DigitalSet & S = shape.pointSet();

  trace.info() << "Skeleton--> "<<S<<std::endl;

  // Display by using two different list to manage OpenGL transparency.
  QApplication application(argc,argv);
  Viewer3D<> viewer;
  viewer.setWindowTitle("homotopicThinning3D");
  viewer.show();

  viewer << SetMode3D( shape_set.className(), "Paving" );
  viewer << CustomColors3D(Color(25,25,255, 255), Color(25,25,255, 255));
  viewer << S ;
  viewer << CustomColors3D(Color(255,25,255, 255), Color(255,25,255, 255));
  viewer << fixedSet;
  viewer << SetMode3D( shape_set.className(), "PavingTransp" );
  viewer << CustomColors3D(Color(250, 0,0, 25), Color(250, 0,0, 5));
  viewer << shape_set;

  viewer<< Viewer3D<>::updateDisplay;
  
  if (vm.count("exportSDP"))
    {
      std::ofstream out;
      out.open(vm["exportSDP"].as<std::string>().c_str());
      for (auto &p : S)
        {
          out << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
    }
  return application.exec();

}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


