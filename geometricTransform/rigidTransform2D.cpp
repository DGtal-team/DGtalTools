 
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

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z2i;
using namespace functors;

/**
@brief Apply rigid transformation on a given image.

@b Usage: rigidTrans2D --input <RawFileName> --output <OutputFileName> --ox 1.0 --oy 1.0 -a 1.2 --tx 1 --ty 0 --m <forward|backward>

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
  --tx FLOAT=0                          X component of translation vector (default 0)
  --ty FLOAT=0                          Y component of translation vector (default 0)

@endcode

@b Example
@code
# transform lena.pgm
./rigidTransform2D -i lena.pgm -o lena_transf.pgm -m backward -a 0.5 --ox 0.5 --oy 0.3 --tx 0.1 --ty 0.1

**/

int main(int argc, char**argv)
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string filename;
  std::string outputFileName;
  std::string model;
  double angle {0};
  double ox {0};
  double oy {0};
  double tx {0};
  double ty {0};

  app.description("Apply rigid transformation on a given image.\n Typical use example:\n \t rigidTrans2D --input <RawFileName> --output <OutputFileName> --ox 1.0 --oy 1.0 -a 1.2 --tx 1 --ty 0 --m <forward|backward>\n");
  app.add_option("-i,--input,1",filename,"Input file.")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output,2",outputFileName,"Output file.")->required();
  app.add_option("-m,--model",model,"Transformation model: backward, forward.")->required();
  app.add_option("-a,--angle",angle,"Rotation angle in radians (default 0)",true);
  app.add_option("--ox",ox,"X coordinate of origin (default 0)",true);
  app.add_option("--oy",oy,"Y coordinate of origin (default 0)",true);
  app.add_option("--tx",tx,"X component of translation vector (default 0)",true);
  app.add_option("--ty",ty,"Y component of translation vector (default 0)",true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

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
