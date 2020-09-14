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
 * @file statisticsEstimators.cpp
 * @ingroup Tools
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr)
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 * @date 2012/06/12
 *
 * Compute statistics between two estimators
 *
 * This file is part of the DGtal tools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <limits>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"

using namespace DGtal;


/**
 @page statisticsEstimators statisticsEstimators
 
 @brief   Computes satistics (L1, L2, Loo) from results of two estimators.

 @b Usage: 	statisticsEstimators --file1 <file1> --column1 <column1> --file2 <file2> --column2 <column2> --output <output>


 @b Allowed @b options @b are : 
 @code
  Positionals:
  1 TEXT:FILE REQUIRED                  File 1.
  2 TEXT:FILE REQUIRED                  File 2.

  Options:
  -h,--help                             Print this help message and exit
  -f,--file1 TEXT:FILE REQUIRED         File 1.
  -F,--file2 TEXT:FILE REQUIRED         File 2.
  -c,--column1 UINT REQUIRED            Column of file 1
  -C,--column2 UINT REQUIRED            Column of file 2
  -o,--output TEXT REQUIRED             Output file
  -m,--monge BOOLEAN=0                  Is from Monge mean computation (optional, default false)
 @endcode

 @b Example: 

 This tool can be used in association to other estimator with for instance the @ref Doc2dLocalEstimators which gives as output a file containing the curvature:

 @code
 ./estimators/2dlocalEstimators --output curvature --shape flower --radius 15 -v 5  --gridstep 1  --estimators 11100 --properties 01
 @endcode

Then you can use this tool as follows:
@code 
 ./estimators/statisticsEstimators --file1 curvature_True_curvature.dat --column1 0 --file2 curvature_II_curvature.dat --column2 0 -o result.dat
@endcode


The resulting file  result.dat should contains:
@verbatim
# h | L1 Mean Error | L2 Mean Error | Loo Mean Error
1 0.0844106 0.0108399 0.519307
@endverbatim
 
 @see
 @ref statisticsEstimators.cpp

 */


/**
 * @brief LoadingStringFromFile Load the current line from an open file, and go to the next line.
 *
 * @param[in,out] file an open file
 * @param[out] value a string who will be override by value of the current line of file
 * @return true if the file is ok, false else
 */
bool LoadingStringFromFile( std::ifstream & file, std::string & value )
{
    if( file.good() )
    {
        std::getline( file, value );
        return true;
    }
    return false;
}

/**
 * @brief split Split a string by a delimiter, and return an array of string contening the splitted string
 *
 * @param[in] s the input string
 * @param[in] delim the delimiter
 * @param[out] elems array filled with splitted strings
 */
void split( const std::string & s, char delim, std::vector< std::string > & elems )
{
    std::stringstream ss( s );
    std::string item;
    while( std::getline( ss, item, delim ))
    {
        elems.push_back( item );
    }
}

/**
 * @brief ComputeStatistics From two file names (string), and two column id, compute statistics (L1, L2, Loo).
 * @param inputdata1 first filename to compare ("myFile1.dat")
 * @param inputdata2 second filename to compare ("myFile2.dat")
 * @param idColumnData1 id of the column to compare from first file
 * @param idColumnData2 id of the column to compare from second file
 * @param isMongeMean hack for some values. Check second file result orientation with first file result. ( file1 : 5, file2 = -5 -> 5 if isMongeMean = true )
 * @param output an iterator who write string
 * @return 1 if all is good, 0 else
 */
int ComputeStatistics ( const std::string & inputdata1,
                        const std::string & inputdata2,
                        const unsigned int & idColumnData1,
                        const unsigned int & idColumnData2,
                        const bool & isMongeMean,
                        std::ofstream & output )
{
    std::ifstream file1( inputdata1.c_str() );
    std::ifstream file2( inputdata2.c_str() );

    double absd1d2;
    double L1 = 0.0;
    double L2 = 0.0;
    double Linf = 0.0;

    std::string s1, s2;
    double v1, v2;
    double h = - std::numeric_limits<double>::max();

    unsigned int nb_elements = 0;
    bool finish = false;
    while(( LoadingStringFromFile( file1, s1 ) && LoadingStringFromFile( file2, s2 )) && !finish )
    {
        while ( s1[ 0 ] == '#' )
        {
            std::size_t  p = s1.find( "# h = " );
            if ( p != std::string::npos )
            {
                h = atof((s1.erase( p, 5 )).c_str());
            }
            if( ! LoadingStringFromFile( file1, s1 ) )
            {
                s1 = "NA";
                finish = true;
            }
        }

        while ( s2[ 0 ] == '#' )
        {
            if( ! LoadingStringFromFile( file2, s2 ) )
            {
                s2 = "NA";
                finish = true;
            }
        }

        if ( s1 == "NA" || s1 == "-nan" || s1 == "-inf" || s1 == "inf" || s1 == "" || s1 == " " )
            continue;
        if ( s2 == "NA" || s2 == "-nan" || s2 == "-inf" || s2 == "inf" || s2 == "" || s2 == " " )
            continue;

        std::vector< std::string > elems1;
        split( s1, ' ', elems1 );
        std::vector< std::string > elems2;
        split( s2, ' ', elems2 );

        if( elems1.size() <= idColumnData1 )
        {
            std::cerr << "Can't found " << idColumnData1 << " column on file1 (" << inputdata1 << "). Is the file/column exist ?" << std::endl;
            continue;
        }
        if( elems2.size() <= idColumnData2 )
        {
            std::cerr << "Can't found " << idColumnData2 << " column on file2 (" << inputdata2 << "). Is the file/column exist ?" << std::endl;
            continue;
        }

        v1 = atof( elems1[ idColumnData1 ].c_str() );
        v2 = atof( elems2[ idColumnData2 ].c_str() );

        if( isMongeMean && (( v1 >= 0.0 ) ^ ( v2 >= 0.0 ))) // hack for Monge. Can be reversed.
        {
            v2 = -v2;
        }

        absd1d2 = std::abs ( v1 - v2 );
        if ( Linf < absd1d2 )
        {
            Linf = absd1d2;
        }
        L1 += absd1d2;
        L2 += absd1d2 * absd1d2;

        ++nb_elements;
    }

    if( h == - std::numeric_limits<double>::max())
    {
        std::cerr << "Can't found h value on file1 (" << inputdata1 << "). Is the file exist ?" << std::endl;
        return 0;
    }

    double meanL1 = L1 / (double)nb_elements;
    double meanL2 = ( sqrt ( L2 )) / (double)nb_elements;
    
    output << h << " "
           << meanL1 << " "
           << meanL2 << " "
           << Linf
           << std::endl;

    return 1;
}

int main( int argc, char** argv )
{
    // parse command line CLI ----------------------------------------------
    CLI::App app;
    std::string filename1;
    std::string filename2;
    unsigned int column1;
    unsigned int column2;
    std::string output_filename;
    bool isMongeMean {false};

    app.description("Computes satistics (L1, L2, Loo) from results of two estimators.\n Typical use example:\n \t statisticsEstimators --file1 <file1> --column1 <column1> --file2 <file2> --column2 <column2> --output <output>\n");
    app.add_option("-f,--file1,1",filename1,"File 1.")->required()->check(CLI::ExistingFile);
    app.add_option("-F,--file2,2",filename2,"File 2.")->required()->check(CLI::ExistingFile);
    app.add_option("--column1,-c", column1, "Column of file 1" )->required();
    app.add_option("--column2,-C", column2, "Column of file 2" )->required();
    app.add_option("--output,-o", output_filename, "Output file")->required();
    app.add_option("--monge,-m", isMongeMean, "Is from Monge mean computation (optional, default false)", true);

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------  
    
    std::ifstream inFileEmptyTest; inFileEmptyTest.open(output_filename.c_str());
    bool isNew = inFileEmptyTest.peek() == std::ifstream::traits_type::eof(); inFileEmptyTest.close();
    std::ofstream file( output_filename.c_str(), std::ofstream::out | std::ofstream::app );

    if( isNew )
    {
        file << "# h | "
             << "L1 Mean Error | "
             << "L2 Mean Error | "
             << "Loo Mean Error"
             << std::endl;
    }

    if ( ComputeStatistics( filename1, filename2, column1, column2, isMongeMean, file ) == 0 )
    {
        file.close();
        return -1;
    }

    file.close();
    return 1;
}
