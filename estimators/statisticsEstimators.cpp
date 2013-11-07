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
#include <vector>
#include <string>
#include <math.h>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

int LoadingStringFromFile ( std::string & parseDataFile, std::string & parseParamFile, const std::string filename )
{
  std::ifstream myfile;
  std::stringstream ss_data;
  std::stringstream ss_param;

  myfile.open ( filename.c_str() );
  if ( myfile.is_open() )
  {
    std::string currentLine;
    while ( myfile.good() )
    {
      getline ( myfile, currentLine );
      if (currentLine.size() > 0 && currentLine[0] == '#')
      {
        ss_param << currentLine << '\n';
      }
      else
      {
        ss_data << currentLine << '\n';
      }
    }
    myfile.close();
    parseDataFile = ss_data.str();
    parseParamFile = ss_param.str();
  }
  else
  {
    //std::cout << "estimatorStatistics@LoadingStringFromFile error : Can't open file " << filename << std::endl;
    return 0; //silent exit
  }
  return 1;
}

int ComputeStatisticsFromString ( const unsigned int idColumnData1, const unsigned int idColumnData2, const std::string & inputdata, const std::string & inputparam )
{

  std::vector<double> data1, data2;

  boost::char_separator<char> sep_lines("\n");
  boost::char_separator<char> sep_column(" ");

  boost::tokenizer< boost::char_separator<char> > tokens_lines(inputdata, sep_lines);
  BOOST_FOREACH (const std::string & line, tokens_lines)
  {
    boost::tokenizer< boost::char_separator<char> > tokens(line, sep_column);
    boost::tokenizer< boost::char_separator<char> >::iterator beg = tokens.begin();
    boost::tokenizer< boost::char_separator<char> >::iterator current = beg;

    boost::tokenizer< boost::char_separator<char> >::iterator itdata1;
    boost::tokenizer< boost::char_separator<char> >::iterator itdata2;

    for ( int i = 0; i < idColumnData1; ++i )
      ++current;
    itdata1 = current;
    current = beg;
    for ( int i = 0; i < idColumnData2; ++i )
      ++current;
    itdata2 = current;

    std::string cstring = ( (std::string)(*itdata1) );
    if ( cstring != "NA" && cstring != "" )
    {
      data1.push_back ( atof ( cstring.c_str() ) );
    }
    else
    {
      return 0;
    }

    cstring = ( (std::string)(*itdata2) );
    if ( cstring != "NA" && cstring != "" )
    {
      data2.push_back ( atof ( cstring.c_str() ) );
    }
    else
    {
      return 0;
    }
  }

  unsigned int sizeVector1 = data1.size();
  unsigned int sizeVector2 = data2.size();
  if ( sizeVector1 == 0 || sizeVector2 == 0 )
  {
    return 0;
  }

  double h = 0.0;
  double radius = 0.0;

  boost::tokenizer< boost::char_separator<char> > tokens_param(inputparam, sep_lines);
  BOOST_FOREACH (const std::string & cline, tokens_param)
  {
    std::string line = cline;
    if ( h != 0.0 && radius != 0.0 )
      break;

    unsigned int pos;

    pos = line.find("# h = ");
    if ( pos != std::string::npos )
    {
      h = atof ( (line.erase ( pos, 5 )).c_str() );
      continue;
    }

    pos = line.find("# computed kernel radius = ");
    if ( pos != std::string::npos)
    {
      radius = atof ( (line.erase ( pos, 26 )).c_str() );
      continue;
    }
  }

  if (h == 0.0 )
  {
    std::cout << "estimatorStatistics@ComputeStatisticsFromString error : h param can't found. h = " << h << std::endl;
    return 0;
  }

  if ( sizeVector1 != sizeVector2 )
  {
    std::cout << "estimatorStatistics@ComputeStatisticsFromString error : data1 & data 2 haven't the same size." << sizeVector1 << " " << sizeVector2 << std::endl;
    return 0;
  }

  double absd1d2 = 0.0;

  double L1 = 0.0;
  double L2 = 0.0;
  double Linf = 0.0;

  for ( int index = 0; index < sizeVector1; ++index )
  {
      absd1d2 = abs ( data1[index] - data2[index] );
      if ( Linf < absd1d2 )
          Linf = absd1d2;
      L1 += absd1d2;
      L2 += absd1d2 * absd1d2;
  }

  std::cout
      << "# h | "
      << "Kernel radius | "
      << "L1 Mean Error | "
      << "L2 Mean Error | "
      << "Loo Mean Error"
      << std::endl;

  double meanL1 = L1 / (double)sizeVector1;
  double meanL2 = ( sqrt ( L2 )) / (double)sizeVector1;
  std::cout
      << h << " "
      << radius << " "
      << meanL1 << " "
      << meanL2 << " "
      << Linf
  << std::endl;

  return 1;
}

int main( int argc, char** argv )
{
  if (argc != 4)
  {
    std::cout << "argc = " << argc << std::endl;
    std::cout << "ERROR. You need to specify a file and two column ID" << std::endl;
  }
  std::string filename = argv[1];
  int column1 = atoi(argv[2]);
  int column2 = atoi(argv[3]);

  std::string inputData;
  std::string inputParam;

  if ( LoadingStringFromFile ( inputData, inputParam, filename ) == 0 )
    return 0;
  if ( ComputeStatisticsFromString ( column1, column2, inputData, inputParam ) == 0 )
    return 0;

  return 0;
}
