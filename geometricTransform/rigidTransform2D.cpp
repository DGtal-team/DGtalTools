 
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
 * @file rigidTransform2D.cpp
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/07/16
 *
 * This file is a part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/images/RigidTransformation2D.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z2i;
using namespace functors;

namespace po = boost::program_options;


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


int main(int argc, char**argv)
{

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input file." )
    ( "output,o", po::value<string>(),"Output filename." )
    ( "model,m", po::value<string>(),"Transformation model: backward, forward." )
    ( "angle,a", po::value<double>(),"Rotation angle in radians." )
    ( "ox", po::value<double>(),"X coordinate of origin." )
    ( "oy", po::value<double>(),"Y coordinate of origin." )
    ( "tx", po::value<double>(),"X component of translation vector." )
    ( "ty", po::value<double>(),"Y component of translation vector." );

  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }

  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Rotate 2D image."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "rigidTrans2D --ox 1.0 --oy 1.0 -a 1.2 --tx 1 --ty 0 --m <forward|backward> --input <RawFileName> --output <VolOutputFileName> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();
  if ( ! ( vm.count ( "model" ) ) ) missingParam ( "--model" );
  std::string model = vm["model"].as<std::string>();
  if ( ! ( vm.count ( "angle" ) ) ) missingParam ( "--angle" );
  double angle =  vm["angle"].as<double>();
  if ( ! ( vm.count ( "ox" ) ) ) missingParam ( "--ox" );
  double ox =  vm["ox"].as<double>();
  if ( ! ( vm.count ( "oy" ) ) ) missingParam ( "--oy" );
  double oy =  vm["oy"].as<double>();
  if ( ! ( vm.count ( "tx" ) ) ) missingParam ( "--tx" );
  double tx =  vm["tx"].as<double>();
  if ( ! ( vm.count ( "ty" ) ) ) missingParam ( "--ty" );
  double ty =  vm["ty"].as<double>();

  typedef ImageSelector<Domain, unsigned char >::Type Image;
  typedef ForwardRigidTransformation2D < Space > ForwardTrans;
  typedef BackwardRigidTransformation2D < Space > BackwardTrans;
  typedef ConstImageAdapter<Image, Domain, BackwardTrans, Image::Value, Identity > MyImageBackwardAdapter;
  typedef DomainRigidTransformation2D < Domain, ForwardTrans > MyTransformedDomain;
  typedef MyTransformedDomain::Bounds Bounds;

  Image image = GenericReader <Image, 2>::import(filename);
  ForwardTrans forwardTrans( RealPoint ( ox, oy ), angle, RealVector( tx, ty ) );
  BackwardTrans backwardTrans( RealPoint ( ox, oy ), angle, RealVector( tx, ty ) );
  MyTransformedDomain domainForwardTrans ( forwardTrans );
  Bounds bounds = domainForwardTrans ( image.domain() );
  Domain transformedDomain ( bounds.first, bounds.second );
  Identity idD;

  if ( model == "forward" )
  {
      Image forwardTransformedImage ( transformedDomain );
      for ( Domain::ConstIterator it = image.domain().begin(); it != image.domain().end(); ++it )
      {
          forwardTransformedImage.setValue ( forwardTrans ( *it ), image ( *it ) );
      }
      GenericWriter<Image, 2, unsigned char, Identity>::exportFile(outputFileName, forwardTransformedImage, idD);
  }
  else
  {
      MyImageBackwardAdapter backwardImageAdapter ( image, transformedDomain , backwardTrans, idD );
      GenericWriter<MyImageBackwardAdapter, 2, unsigned char, Identity>::exportFile(outputFileName, backwardImageAdapter, idD);
  }
  return 0;
}
