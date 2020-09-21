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
 * @file longvol2vol.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/05/01
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/readers/LongvolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/base/BasicFunctors.h>

#include "CLI11.hpp"

using namespace DGtal;
using namespace Z3i;


/**
 @page longvol2vol longvol2vol
 @brief  Converts a longvol (long int) to a vol file (unsigned char).

@b Usage: longvol2vol [input] [output]

@b Allowed @b options @b are:

@code
Positionals:
  1 TEXT:FILE REQUIRED                  Input longvol filename ( .longvol)
  2 TEXT                                Output vol filename.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input longvol filename ( .longvol)
  -o,--output TEXT                      Output vol filename.
  -m,--mode UINT:{0,1,2}                Conversion mode:
                                        	 0 = cast (default)
                                        	 1 = Linear Scaling
                                        	 2 = Grayscale cycle (32 steps, except 0 values).
@endcode

@b Example:
@code
  $ longvol2vol -i ${DGtal}/tests/samples/test.longvol -o out.vol 
@endcode

@see longvol2vol.cpp

*/

/**
 * Cycle scaling functor
 */
struct CycleFunctor
{
  template<typename TInput>
  inline
  unsigned char operator()(const TInput& aInput) const
  {
     if (aInput == 0)
       return 0;
     else
       return static_cast<unsigned char>(1+NumberTraits<TInput>::castToInt64_t(aInput) % 32);
  }
};


/**
 * Linear scaling functor
 */
struct LinearFunctor
{
  LinearFunctor(const DGtal::uint64_t vmax, const DGtal::uint64_t vmin):
  myMin(vmin), myMax(vmax) {}
  
  inline unsigned char operator()(const DGtal::int64_t& aInput) const
  {
    return static_cast<unsigned char>( (NumberTraits<DGtal::uint64_t>::castToDouble(aInput - myMin)*255.0 /
                                        (double) (myMax-myMin)));
  }
  
  DGtal::uint64_t myMin,myMax;
};


int main(int argc, char**argv)
{
  

   // parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.vol"};
   unsigned int mode {0};
   
   app.description("Converts a longvol (long int) to a vol file (unsigned char). \n Basic example:\n\t longvol2vol --input <LongvolFileName> --o <VolOutputFileName> ");
   app.add_option("-i,--input,1", inputFileName, "Input longvol filename ( .longvol)" )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2",outputFileName,"Output vol filename." );
   app.add_option("-m,--mode", mode, "Conversion mode:\n\t 0 = cast (default)\n\t 1 = Linear Scaling\n\t 2 = Grayscale cycle (32 steps, except 0 values).")
     -> check(CLI::IsMember({0, 1, 2}));

  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  //Main program
  typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t>  MyImageC;
  MyImageC  imageC = LongvolReader< MyImageC >::importLongvol ( inputFileName );
  bool res = false;
  
  if (mode == 0)
    res =  VolWriter< MyImageC, functors::Cast<unsigned char> >::exportVol(outputFileName, imageC);
  else
    if (mode == 1)
    {
      DGtal::uint64_t vmax, vmin;
      vmax = *std::max_element(imageC.begin(), imageC.end());
      vmin = *std::min_element(imageC.begin(), imageC.end());
      trace.info() << "Max value = "<< vmax<<std::endl;
      trace.info() << "Min value = "<< vmin<<std::endl;
      res  = VolWriter<MyImageC , LinearFunctor>::exportVol(outputFileName, imageC, true, LinearFunctor(vmax,vmin));
    }
  if (mode == 2)
    res = VolWriter<MyImageC, CycleFunctor>::exportVol(outputFileName, imageC);
  
  if (res)
    return EXIT_SUCCESS;
  
  else
  {
    trace.error()<<"Error while exporting the Vol.";
    trace.info()<<std::endl;
    return EXIT_FAILURE;
  }
}
