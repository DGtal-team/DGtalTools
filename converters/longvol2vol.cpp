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
 * @file vol2raw.cpp
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
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

/**
 @page longvol2vol longvol2vol
 @brief  Converts a longvol (long int) to a vol file (unsigned char).

@b Usage: longvol2vol [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]          display this message.
  -i [ --input ] arg     Input longvol filename.
  -o [ --output ] arg    Output vole filename.
  -m [ --mode ] arg (=0) Conversion mode:
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
  ( "input,i", po::value<std::string>(), "Input longvol filename." )
  ( "output,o", po::value<std::string>(),"Output vole filename." )
  ( "mode,m", po::value<unsigned int>()->default_value(0),"Conversion mode:\n\t 0 = cast (default)\n\t 1 = Linear Scaling\n\t 2 = Grayscale cycle (32 steps, except 0 values)." )
  ;
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ) ||argc<=1 )
  {
    trace.info() << "Converts a longvol (long int) to a vol file (unsigned char)."<<std::endl
    << std::endl << "Basic usage: "<<std::endl
    << "\tlongvol2vol --input <LongvolFileName> --o <VolOutputFileName> "<<std::endl
    << general_opt << "\n";
    return 0;
  }
  
  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();
  unsigned int mode = vm["mode"].as<unsigned int>();
  
  
  //Main program
  typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t>  MyImageC;
  MyImageC  imageC = LongvolReader< MyImageC >::importLongvol ( filename );
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
    return 0;
  else
  {
    trace.error()<<"Error while exporting the Vol.";
    trace.info()<<std::endl;
    return 1;
  }
}
