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

#include "CLI11.hpp"

#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;


///////////////////////////////////////////////////////////////////////////////

/**
 @page homotopicThinning3D homotopicThinning3D
 
 @brief Applies an homotopic thinning of a 3d image file (vol,longvol,pgm3d...) with 3D viewer.

 @b Usage: homotopicThinning3d [options] --input <3dImageFileName>  {vol,longvol,pgm3d...}


 @b Allowed @b options @b are : 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  Input volumetric file (.vol, .pgm3d or p3d)

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input volumetric file (.vol, .pgm3d or p3d)
   -m,--min FLOAT=0                      Minimum (excluded) value for threshold.
   -M,--max FLOAT=255                    Maximum (included) value for threshold.
   -e,--exportSDP TEXT                   Export the resulting set of points in a simple (sequence of discrete point (sdp)).
   --fixedPoints INT ...                 defines the coordinates of points which should not be removed.
   -s,--fixedPointSDP TEXT:FILE          use fixed points from a file.
   
 
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
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string exportSDPName;
  std::string fixedPointSDPName;
  double min {0};
  double max {255};
  
  std::vector<int> vectFixedPoints;
  
  app.description("Applies an homotopic thinning of a 3d image file (vol,longvol,pgm3d...) with 3D viewer. \n Basic usage: \t homotopicThinning3d [options] --input <3dImageFileName>  {vol,longvol,pgm3d...}  \n Usage by forcing point to be left by the thinning: \n homotopicThinning3D --input ${DGtal}/examples/samples/Al.100.vol  --fixedPoints 56 35 5  56 61 5  57 91 38  58 8 38  45 50 97 \n");
  app.add_option("-i,--input,1", inputFileName, "Input volumetric file (.vol, .pgm3d or p3d)" )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("--min,-m", min, "Minimum (excluded) value for threshold.", true);
  app.add_option("--max,-M", max, "Maximum (included) value for threshold.", true);
  app.add_option("--exportSDP,-e",exportSDPName, "Export the resulting set of points in a simple (sequence of discrete point (sdp)).");
  app.add_option("--fixedPoints",vectFixedPoints, "defines the coordinates of points which should not be removed.");
  app.add_option("--fixedPointSDP,-s",fixedPointSDPName,  "use fixed points from a file.")
   ->check(CLI::ExistingFile);
  
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  typedef ImageSelector < Z3i::Domain, unsigned char>::Type Image;
  Image image = GenericReader<Image>::import ( inputFileName );

  trace.beginBlock("DT Computation");
  typedef functors::IntervalForegroundPredicate<Image> Predicate;
  Predicate aPredicate(image, min, max );

  const Z3i::L2Metric aMetric{};
  DistanceTransformation<Z3i::Space, Predicate , Z3i::L2Metric> dt(image.domain(), aPredicate, aMetric );
  trace.endBlock();
  trace.info() <<image<<std::endl;

  // Domain creation from two bounding points.
  trace.beginBlock("Constructing Set");
  DigitalSet shape_set( image.domain() );
  DigitalSet fixedSet( image.domain() );

  // Get the optional fixed points
  if(vectFixedPoints.size() > 0){
    if(vectFixedPoints.size()%3==0){
      for( unsigned int i=0; i < vectFixedPoints.size()-2; i=i+3)
        {
          Z3i::Point pt(vectFixedPoints.at(i), vectFixedPoints.at(i+1), vectFixedPoints.at(i+2));
          fixedSet.insertNew(pt);
        }
    }
    else
    {
      trace.error()<< " The coordinates should be 3d coordinates, ignoring fixedPoints option." << std::endl;
    }
  }
  
  if(fixedPointSDPName != "")
    {
      std::vector<Z3i::Point> vPt = PointListReader<Z3i::Point>::getPointsFromFile(fixedPointSDPName);
      for( auto &p: vPt)
        {
          fixedSet.insert(p);
        }
    }
  
  SetFromImage<DigitalSet>::append<Image>(shape_set, image, min, max );
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
  
  if (exportSDPName != "")
    {
      std::ofstream out;
      out.open(exportSDPName.c_str());
      for (auto &p : S)
        {
          out << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
    }
  return application.exec();

}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


