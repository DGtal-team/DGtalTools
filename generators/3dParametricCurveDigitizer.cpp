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
 * @file 3dParametricCurveDigitizer.cpp
 * @ingroup Tools
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr)
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2018/08/16
 *
 * DGtal 3D parametric curve digitier.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <string>
#include <iterator>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#include "DGtal/geometry/curves/parametric/EllipticHelix.h"
#include "DGtal/geometry/curves/parametric/Knot_3_1.h"
#include "DGtal/geometry/curves/parametric/Knot_3_2.h"
#include "DGtal/geometry/curves/parametric/Knot_4_1.h"
#include "DGtal/geometry/curves/parametric/Knot_4_3.h"
#include "DGtal/geometry/curves/parametric/Knot_5_1.h"
#include "DGtal/geometry/curves/parametric/Knot_5_2.h"
#include "DGtal/geometry/curves/parametric/Knot_6_2.h"
#include "DGtal/geometry/curves/parametric/Knot_7_4.h"
#include "DGtal/geometry/curves/parametric/NaiveParametricCurveDigitizer3D.h"
#include "DGtal/geometry/curves/parametric/DecoratorParametricCurveTransformation.h"
#include "DGtal/images/RigidTransformation3D.h"

using namespace DGtal;


/**
 @page 3dParametricCurveDigitizer 3dParametricCurveDigitizer
 @brief  Digitizes 3D parametric curves using DGtal library.
 


 @b Usage:  3dParametricCurveDigitizer [options] --curve <curve> --param1 <double> --param2 <double> --param3 <double> --tstart <double> --tend <double> --step <double> --output <basename>

 @b Allowed @b options @b are:

 @code
  -h,--help                             Print this help message and exit
  --param1 FLOAT=1                      a radius or a scaling factor (default 0)
  --param2 FLOAT=1                      a radius or a scaling factor (default 0)
  --param3 FLOAT=1                      a radius or a scaling factor (default 0)
  --tstart FLOAT                        start time
  --tend FLOAT                          end time
  -s,--step FLOAT                       step
  -k,--knext UINT=5                     K_NEXT value (default 5)
  -l,--list                             List all available shapes
  -c,--curve TEXT                       Shape name
  -a,--angle FLOAT=0                    Rotation angle in radians(default 0)
  --ox FLOAT=0                          X coordinate of origin (default 0)
  --oy FLOAT=0                          Y coordinate of origin (default 0)
  --oz FLOAT=0                          Z coordinate of origin (default 0)
  --ax FLOAT=1                          X component of rotation axis (default 1)
  --ay FLOAT=0                          Y component of rotation axis (default 0)
  --az FLOAT=0                          Z component of rotation axis (default 0)
  -o,--output TEXT                      Basename of the output file

 @endcode
 You can list the potential curves:
 @code
 $ 3dParametricCurveDigitizer --list
 3D Parametric curves:
	EllipticHelix Helix with two axes and a distance between consecutive turns. Its period is 2 * Pi.
		Required parameter(s): --param1, --param2, --param3
	Knot_3_1 Parametric knot 3_1 i.e., a trefoli polynomial knot
		Required parameter(s): --param1, --param2, --param3
	Knot_3_2 Parametric knot 3_2 i.e., a trefoli knot having a period of 2 * Pi
		Required parameter(s): --param1, --param2, --param3
	Knot_4_1 Parametric knot 4_1 (polynomial knot)
		Required parameter(s): --param1, --param2, --param3
	Knot_4_3 Parametric knot 4_3 having a period 2 * Pi.
		Required parameter(s): --param1, --param2, --param3
	Knot_5_1 Parametric knot 5_1 (polynomial knot)
		Required parameter(s): --param1, --param2, --param3
	Knot_5_2 Parametric knot 5_2 having a period 2 * Pi
		Required parameter(s): --param1, --param2, --param3
	Knot_6_2 Parametric knot 6_2 (polynomial knot)
		Required parameter(s): --param1, --param2, --param3
	Knot_7_4 Parametric knot 7_4 (polynomial knot)
		Required parameter(s): --param1, --param2, --param3
 @endcode

 @b Example:
 @code
 # generate a trefoli knot:
 3dParametricCurveDigitizer --curve Knot_3_1 --param1 10 --param2 10 --param3 10 --tstart -2.2 --tend 2.2 --step 0.001 --output knot_3_1
 @endcode


 @see
 @ref 3dParametricCurveDigitizer.cpp
 @ref 3dParametricCurveDigitizer
 */

/**
 * Global vectors to describe the available curves.
 */
std::vector<std::string> curves3D;
std::vector<std::string> curvesDesc;

/**
 * Create the static list of shapes.
 *
 */
void createList()
{
  curves3D.push_back("EllipticHelix");
  curvesDesc.push_back("A helix with two axes.");

  curves3D.push_back("Knot_3_1");
  curvesDesc.push_back("A parametric polynomial knot.");

  curves3D.push_back("Knot_3_2");
  curvesDesc.push_back("A parametric knot having period 2 * Pi.");

  curves3D.push_back("Knot_4_1");
  curvesDesc.push_back("A parametric polynomial knot.");

  curves3D.push_back("Knot_4_3");
  curvesDesc.push_back("A parametric knot having period 2 * Pi.");

  curves3D.push_back("Knot_5_1");
  curvesDesc.push_back("A parametric polynomial knot.");

  curves3D.push_back("Knot_5_2");
  curvesDesc.push_back("A parametric knot having period 2 * Pi.");

  curves3D.push_back("Knot_6_2");
  curvesDesc.push_back("A parametric polynomial knot.");

  curves3D.push_back("Knot_7_4");
  curvesDesc.push_back("A parametric polynomial knot.");
}

/**
 * Display the shape list with parameters.
 *
 */
void displayList()
{
  trace.emphase()<<"3D Curves:"<<std::endl;
  for ( unsigned int i = 0; i < curves3D.size(); ++i )
    trace.info()<<"\t"<<curves3D[i]<<"\t"  <<curvesDesc[i]<<std::endl;
}


/**
 * Check if a given curve is available. If not, we exist with an error.
 * If it is, we return the corresponding index in the global vectors.
 *
 * @param curveName name of the curve to search.
 *
 * @return index of the curve in the curve vectors.
 */
unsigned int checkAndRetrunIndex(const std::string &curveName)
{
  unsigned int pos = 0;

  while ((pos < curves3D.size()) && (curves3D[pos] != curveName))
    pos++;

  if ( pos == curves3D.size ( ) )
  {
    trace.error ( ) << "The specified shape has not found.";
    trace.info ( ) << std::endl;
    exit ( 1 );
  }

  return pos;
}


/**
 * Functor to export a given shape into an image file
 * (pgm,raw,pdf,svg,...) and to extract its signature.
 *
 * @tparam Set type of the input Set
 * @tparam Image type of the input Image.
 */
template <typename Iterator>
struct Exporter
{

  /**
   * Export a given Set into an image file.
   *
   * @param begin - begin iterator.
   * @param end - end iterator.
   * @param outputName output file name.
   *
   */
  static
  void save ( Iterator begin, Iterator end, const std::string outputName )
  {
    std::ofstream outfile;
    outfile.open ( outputName + ".dat", std::ofstream::out );
    for ( auto it = begin; it != end; ++it )
      outfile << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;
  }
};

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info()<<std::endl;
  exit(1);
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  std::string curveName;
  std::string outputName;
  double param1 {1};
  double param2 {1};
  double param3 {1};
  double tstart;
  double tend;
  double step;
  unsigned int knext {5};
  double angle;
  double ox {0}, oy {0}, oz {0};
  double ax {1}, ay {0}, az {0}; 

  app.description("Digitizes 3D parametric curves using DGtal library.\n Typical use example:\n \t 3dParametricCurveDigitizer [options] --curve <curve> --param1 <double> --param2 <double> --param3 <double> --tstart <double> --tend <double> --step <double> --output <basename>\n");
  app.add_option("--param1",param1,"a radius or a scaling factor (default 0)",true);
  app.add_option("--param2",param2,"a radius or a scaling factor (default 0)",true);
  app.add_option("--param3",param3,"a radius or a scaling factor (default 0)",true);
  auto tstartOpt = app.add_option("--tstart",tstart,"start time");
  auto tendOpt = app.add_option("--tend",tend,"end time");
  auto stepOpt = app.add_option("--step, -s",step,"step");
  app.add_option("--knext, -k",knext,"K_NEXT value (default 5)",true);
  auto listOpt = app.add_flag("--list,-l","List all available shapes");
  auto curveNameOpt = app.add_option("--curve,-c",curveName,"Shape name");
  app.add_option("--angle,-a",angle,"Rotation angle in radians(default 0)",true);
  app.add_option("--ox",ox,"X coordinate of origin (default 0)",true);
  app.add_option("--oy",oy,"Y coordinate of origin (default 0)",true);
  app.add_option("--oz",oz,"Z coordinate of origin (default 0)",true);
  app.add_option("--ax",ax,"X component of rotation axis (default 1)",true);
  app.add_option("--ay",ay,"Y component of rotation axis (default 0)",true);
  app.add_option("--az",az,"Z component of rotation axis (default 0)",true);
  auto outputNameOpt = app.add_option("--output,-o",outputName,"Basename of the output file");

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  //List creation
  createList();

  if ( listOpt->count() > 0 )
  {
    displayList();
    return 0;
  }

  if ( curveNameOpt->count() == 0) missingParam("--curve");
  if ( outputNameOpt->count() == 0) missingParam("--output");
  if ( tstartOpt->count() == 0) missingParam("--tstart");
  if ( tendOpt->count() == 0) missingParam("--tend");
  if ( stepOpt->count() == 0) missingParam("--step");

  typedef functors::ForwardRigidTransformation3D < Z3i::Space, Z3i::RealPoint, Z3i::RealPoint, functors::Identity > ForwardTrans;
  ForwardTrans trans ( Z3i::RealPoint ( ox, oy, oz ), Z3i::RealPoint ( ax, ay, az ), angle, Z3i::RealVector ( 0, 0, 0 ) );
  //We check that the shape is known
  unsigned int id = checkAndRetrunIndex ( curveName );

  if ( id == 0 )
  {
    typedef EllipticHelix < Z3i::Space > MyHelix;
    typedef DecoratorParametricCurveTransformation < MyHelix, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef typename NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyHelix helix ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( helix, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 1 )
  {
    typedef Knot_3_1 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 2 )
  {
    typedef Knot_3_2 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 3 )
  {
    typedef Knot_4_1 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 4 )
  {
    typedef Knot_4_3 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 5 )
  {
    typedef Knot_5_1 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 6 )
  {
    typedef Knot_5_2 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 7 )
  {
    typedef Knot_6_2 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  if ( id == 8 )
  {
    typedef Knot_7_4 < Z3i::Space > MyKnot;
    typedef DecoratorParametricCurveTransformation < MyKnot, ForwardTrans > MyRotatedCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef NaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

    MyDigitalCurve digitalCurve;
    MyKnot knot ( param1, param2, param3 );
    MyRotatedCurve rotCurve ( knot, trans );
    Digitizer digitize;
    digitize.init ( tstart, tend, step );
    digitize.setKNext ( knext );
    digitize.attach ( &rotCurve );
    digitize.digitize ( std::back_insert_iterator < MyDigitalCurve> ( digitalCurve ) );
    Exporter < MyDigitalCurve::const_iterator >::save ( digitalCurve.begin ( ), digitalCurve.end ( ), outputName );
  }

  return 0;
}
