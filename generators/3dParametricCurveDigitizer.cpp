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

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <icc34.h>
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
#include "DGtal/geometry/curves/parametric/UglyNaiveParametricCurveDigitizer3D.h"
#include "DGtal/geometry/curves/parametric/DecoratorParametricCurveTransformation.h"
#include "DGtal/images/RigidTransformation3D.h"

using namespace DGtal;


/**
 @page shapeGenerator shapeGenerator
 @brief  Generates shapes using DGtal library.
 


 @b Usage:  shapeGenerator [options] --shape <shapeName> --output <outputBasename>

 @b Allowed @b options @b are:

 @code
  -h  [ --help ]                   Display this message
  -p1 [ --param1]                  First parameter e.g. a radius or a scaling factor
  -p2 [ --param2]                  Second parameter e.g. a radius or a scaling factor
  -p3 [ --param3]                  Third parameter e.g. a distance between consecutive turns or a scaling factor
  -ts [ --tstart]                  Start time
  -te [ --tend]                    End time
  -s  [ --step]                    Step
  -k  [ --knext]                   K_NEXT value
  -l  [ --list ]                   List all available curves
  -c  [ --curve]                   A curve to digitize
  -a [--angle],                    Rotation angle in radians
  -ox                              X coordinate of origin
  -oy                              Y coordinate of origin
  -oz                              Z coordinate of origin
  -ax                              X component of rotation axis
  -ay                              Y component of rotation axis
  -az                              Z component of rotation axis
  -o  [ --output ] arg             Basename of the output file

 @endcode
 You can list the potential curves:
 @code
 $ contourGenerator --list
 3D Parametric curves:
	EllipticHelix Helix with two axes and a distance between consecutive turns. Its period is 2 * Pi.
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_3_1 Parametric knot 3_1 i.e., a trefoli polynomial knot
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_3_2 Parametric knot 3_2 i.e., a trefoli knot having a period of 2 * Pi
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_4_1 Parametric knot 4_1 (polynomial knot)
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_4_3 Parametric knot 4_3 having a period 2 * Pi.
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_5_1 Parametric knot 5_1 (polynomial knot)
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_5_2 Parametric knot 5_2 having a period 2 * Pi
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_6_2 Parametric knot 6_2 (polynomial knot)
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
	Knot_7_4 Parametric knot 7_4 (polynomial knot)
		Required parameter(s): --param1 [-p1], --param2 [p2], --param3 [-p3]
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
void missingParam ( std::string param )
{
  trace.error ( ) <<" Parameter: "<< param <<" is required..";
  trace.info ( ) << std::endl;
  exit ( 1 );
}

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("param1,p1", po::value<double>()->default_value(1), "a radius or a scaling factor")
  ("param2,p2", po::value<double>()->default_value(1), "a radius or a scaling factor")
  ("param3,p3", po::value<double>()->default_value(1), "a radius or a scaling factor")
  ("tstart,ts", po::value<double>(), "start time")
  ("tend,te", po::value<double>(), "end time")
  ("step,s", po::value<double>(), "step")
  ("knext,k",  po::value<unsigned int>()->default_value(5), "K_NEXT value" )
  ("list,l",  "List all available shapes")
  ("curve,c", po::value<std::string>(), "Shape name")
  ( "angle,a", po::value<double>()->default_value(0),"Rotation angle in radians." )
  ( "ox", po::value<double>()->default_value(0),"X coordinate of origin." )
  ( "oy", po::value<double>()->default_value(0),"Y coordinate of origin." )
  ( "oz", po::value<double>()->default_value(0),"Z coordinate of origin." )
  ( "ax", po::value<double>()->default_value(1),"X component of rotation axis." )
  ( "ay", po::value<double>()->default_value(0),"Y component of rotation axis." )
  ( "az", po::value<double>()->default_value(0),"Z component of rotation axis." )
  ("output,o", po::value<std::string>(), "Basename of the output file");

  bool parseOK = true;
  po::variables_map vm;
  try
  {
    po::store ( po::parse_command_line ( argc, argv, general_opt ), vm );
  }
  catch ( const std::exception& ex )
  {
    parseOK = false;
    trace.info ( ) << "Error checking program options: "<< ex.what ( ) << std::endl;
  }

  po::notify ( vm );
  if ( ! parseOK || vm.count("help") || argc <= 1 )
  {
    trace.info()<< "Digitizes curves using DGtal library" <<std::endl << "Basic usage: "<<std::endl
    << "\t3dParametricCurveDigitize [options] --curve <curveName> --tstart <double> --tend <double> --step <double> --knext <unsigned int> --output <outputBasename>"<<std::endl
    << general_opt << "\n";
    return 0;
  }

  //List creation
  createList();

  if ( vm.count ( "list" ) )
  {
    displayList();
    return 0;
  }

  //Parse options
  if ( !( vm.count ( "curve" ) ) ) missingParam ( "--curve" );
  std::string curveName = vm["curve"].as < std::string > ( );

  if ( !( vm.count ( "param1" ) ) ) missingParam ( "--param1" );
  double param1 = vm["param1"].as < double > ( );

  if ( !( vm.count ( "param2" ) ) ) missingParam ( "--param2" );
  double param2 = vm["param2"].as < double > ( );

  if ( !( vm.count ( "param3" ) ) ) missingParam ( "--param3" );
  double param3 = vm["param3"].as < double > ( );

  if ( !( vm.count ( "tstart" ) ) ) missingParam ( "--tstart" );
  double tstart = vm["tstart"].as < double > ( );

  if ( !( vm.count ( "tend" ) ) ) missingParam ( "--tend" );
  double tend = vm["tend"].as < double > ( );

  if ( !( vm.count ( "step" ) ) ) missingParam ( "--step" );
  double step = vm["step"].as < double > ( );

  unsigned int knext = vm["knext"].as < unsigned int > ( );

  if ( !( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputName = vm["output"].as < std::string > ( );

  double angle =  vm["angle"].as<double>();
  double ox =  vm["ox"].as<double>();
  double oy =  vm["oy"].as<double>();
  double oz =  vm["oz"].as<double>();
  double ax =  vm["ax"].as<double>();
  double ay =  vm["ay"].as<double>();
  double az =  vm["az"].as<double>();


  typedef functors::ForwardRigidTransformation3D < Z3i::Space, functors::Identity, Z3i::RealPoint, Z3i::RealPoint > ForwardTrans;
  ForwardTrans trans ( Z3i::RealPoint ( ox, oy, oz ), Z3i::RealPoint ( ax, ay, az ), angle, Z3i::RealVector ( 0, 0, 0 ) );
  //We check that the shape is known
  unsigned int id = checkAndRetrunIndex ( curveName );

  if ( id == 0 )
  {
    typedef EllipticHelix < Z3i::Space > MyHelix;
    typedef DecoratorParametricCurveTransformation < MyHelix, ForwardTrans > MyRotatedCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef typename UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >  Digitizer;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::DigitalCurve MyDigitalCurve;
    typedef UglyNaiveParametricCurveDigitizer3D < MyRotatedCurve >::MetaData MyMetaData;

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
