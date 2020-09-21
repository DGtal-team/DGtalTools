 
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
 * @file rigidTransforma3D.cpp
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
#include <DGtal/images/RigidTransformation3D.h>

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;
using namespace functors;

/**
@brief Apply rigid transformation on a given volumic image.

@b Usage : rigidTrans3D --input <RawFileName> --output <VolOutputFileName> --ox 1.0 --oy 1.0 --oz 1 -a 1.2 --ax 1 --ay 1 --az 0 --tx 1 --ty 0 --tz 0 --m <forward|backward>  

@b Allowed @b options @b are :

@code
Positionals:
  1 TEXT:FILE REQUIRED                  Input file.
  2 TEXT REQUIRED                       Output file.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input file.
  -o,--output TEXT REQUIRED             Output file.
  -m,--model TEXT REQUIRED              Transformation model: backward, forward.
  -a,--angle FLOAT=0                    Rotation angle in radians (default 0)
  --ox FLOAT=0                          X coordinate of origin (default 0)
  --oy FLOAT=0                          Y coordinate of origin (default 0)
  --oz FLOAT=0                          Z coordinate of origin (default 0)
  --ax FLOAT=1                          X component of rotation axis (default 1)
  --ay FLOAT=0                          Y component of rotation axis (default 0)
  --az FLOAT=0                          Z component of rotation axis (default 0)
  --tx FLOAT=0                          X component of translation vector (default 0)
  --ty FLOAT=0                          Y component of translation vector (default 0)
  --tz FLOAT=0                          Y component of translation vector (default 0)

@endcode

@b Example

*/

int main(int argc, char**argv)
{

  // parse command line CLI ----------------------------------------------
  CLI::App app;
  std::string filename;
  std::string outputFileName;
  std::string model;
  double angle {0};
  double ox {0}, oy {0}, oz {0};
  double ax {1}, ay {0}, az {0};
  double tx {0}, ty {0}, tz {0};

  app.description("Apply rigid transformation on a given volumic image.\n Typical use example:\n \t rigidTrans3D --input <RawFileName> --output <VolOutputFileName> --ox 1.0 --oy 1.0 --oz 1 -a 1.2 --ax 1 --ay 1 --az 0 --tx 1 --ty 0 --tz 0 --m <forward|backward>\n");
  app.add_option("-i,--input,1",filename,"Input file.")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output,2",outputFileName,"Output file.")->required();
  app.add_option("-m,--model",model,"Transformation model: backward, forward.")->required();
  app.add_option("-a,--angle",angle,"Rotation angle in radians (default 0)",true);
  app.add_option("--ox",ox,"X coordinate of origin (default 0)",true);
  app.add_option("--oy",oy,"Y coordinate of origin (default 0)",true);
  app.add_option("--oz",oz,"Z coordinate of origin (default 0)",true);
  app.add_option("--ax",ax,"X component of rotation axis (default 1)",true);
  app.add_option("--ay",ay,"Y component of rotation axis (default 0)",true);
  app.add_option("--az",az,"Z component of rotation axis (default 0)",true);
  app.add_option("--tx",tx,"X component of translation vector (default 0)",true);
  app.add_option("--ty",ty,"Y component of translation vector (default 0)",true);
  app.add_option("--tz",tz,"Y component of translation vector (default 0)",true);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  typedef ImageSelector<Domain, unsigned char >::Type Image;
  typedef ForwardRigidTransformation3D < Space > ForwardTrans;
  typedef BackwardRigidTransformation3D < Space > BackwardTrans;
  typedef ConstImageAdapter<Image, Domain, BackwardTrans, Image::Value, Identity > MyImageBackwardAdapter;
  typedef DomainRigidTransformation3D < Domain, ForwardTrans > MyTransformedDomain;
  typedef MyTransformedDomain::Bounds Bounds;

  Image image = GenericReader <Image, 3>::import(filename);
  ForwardTrans forwardTrans( RealPoint ( ox, oy, oz ), RealPoint ( ax, ay, az ), angle, RealVector( tx, ty, tz ) );
  BackwardTrans backwardTrans( RealPoint ( ox, oy, oz ), RealPoint ( ax, ay, az ), angle, RealVector( tx, ty, tz ) );
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
      GenericWriter<Image, 3, unsigned char, Identity>::exportFile(outputFileName, forwardTransformedImage," ", idD);
  }
  else
  {
      MyImageBackwardAdapter backwardImageAdapter ( image, transformedDomain , backwardTrans, idD );
      GenericWriter<MyImageBackwardAdapter, 3, unsigned char, Identity>::exportFile(outputFileName, backwardImageAdapter, " ", idD);
  }
  return 0;
}
