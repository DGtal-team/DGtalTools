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
 * @file qglViewer.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <QImageReader>
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/Color.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
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
    ( "input,i", po::value<std::string>(), "Input vol file." );
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
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  
  
  typedef ImageSelector < Z3i::Domain, unsigned char>::Type Image;
  Image image = GenericReader<Image>::import ( filename );

  trace.beginBlock("DT Computation");
  typedef SimpleThresholdForegroundPredicate<Image> Predicate;
  Predicate aPredicate(image,0);

  DistanceTransformation<Z3i::Space, Predicate , Z3i::L2Metric> dt(image.domain(),aPredicate, Z3i::L2Metric() );
  trace.endBlock();
  trace.info() <<image<<std::endl;

  // Domain cretation from two bounding points.
  Point c( 0, 0, 0 );
  Point p1( -50, -50, -50 );
  Point p2( 50, 50, 50 );
  Domain domain( p1, p2 );
  
  trace.beginBlock("Constructing Set");
  DigitalSet shape_set( domain );
  SetFromImage<DigitalSet>::append<Image>(shape_set, image,
                                          0, 255);
  trace.info() << shape_set<<std::endl;
  trace.endBlock();

  trace.beginBlock("Computing skeleton");
  Object26_6 shape( dt26_6, shape_set );
  int nb_simple=0; 
  int layer = 1;
  std::queue<DigitalSet::Iterator> Q;
  do 
    {
      trace.info() << "Layer: "<< layer << std::endl;
      int nb=0;
      DigitalSet & S = shape.pointSet();
 
      for ( DigitalSet::Iterator it = S.begin(); it != S.end(); ++it )
        {
	  trace.progressBar((double)nb, (double)S.size()); 
          nb++;
	  if (dt( *it ) <= layer)
	    {
	      if ( shape.isSimple( *it ) )
		Q.push( it );
	    }
	}
      nb_simple = 0;
      while ( ! Q.empty() )
        {
          DigitalSet::Iterator it = Q.front();
          Q.pop();
          if ( shape.isSimple( *it ) )
            {
              S.erase( *it );
              ++nb_simple;
            }
        }
      trace.info() << "Nb simple points : "<<nb_simple<<std::endl;
      ++layer;
     }
  while ( nb_simple != 0 );
  trace.endBlock();
  
  DigitalSet & S = shape.pointSet();

  trace.info() << "Skeleton--> "<<S<<std::endl;

  // Display by using two different list to manage OpenGL transparency.
  QApplication application(argc,argv);
  Viewer3D<> viewer;
  viewer.setWindowTitle("simpleExample3DViewer");
  viewer.show();  
  
  viewer << SetMode3D( shape_set.className(), "Paving" );
  viewer << CustomColors3D(Color(25,25,255, 255), Color(25,25,255, 255));
  viewer << S ; 

  viewer << SetMode3D( shape_set.className(), "PavingTransp" );
  viewer << CustomColors3D(Color(250, 0,0, 25), Color(250, 0,0, 5));
  viewer << shape_set;

  viewer<< Viewer3D<>::updateDisplay;
   
  return application.exec();

}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


