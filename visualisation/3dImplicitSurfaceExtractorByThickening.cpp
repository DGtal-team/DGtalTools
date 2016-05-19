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
 * @file 3dImplicitSurfaceExtractorByThickening.cpp
 * @ingroup Visualisation
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/08/28
 *
 * A tool to visualize 3D implicit surface by thickening the zero
 * subset and then shrink it with cubical complex collapses.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <map>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/base/CountedPtr.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/CubicalComplex.h"
#include "DGtal/topology/CubicalComplexFunctions.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;



/**
 @page  Doc3dImplicitSurfaceExtractorByThickening 3dImplicitSurfaceExtractorByThickening
 
 @brief Computes the zero level set of the given polynomial.

 @b Usage:  3dImplicitSurfaceExtractorByThickening [options] input

 @b Allowed @b options @b are :
 
 @code
  -h [ --help ]                         display this message
  -p [ --polynomial ] arg               the implicit polynomial whose 
                                        zero-level defines the shape of 
                                        interest.
  -a [ --minAABB ] arg (=-10)           the min value of the AABB bounding box 
                                        (domain)
  -A [ --maxAABB ] arg (=10)            the max value of the AABB bounding box 
                                        (domain)
  -g [ --gridstep ] arg (=1)            the gridstep that defines the 
                                        digitization (often called h). 
  -t [ --thickness ] arg (=0.01)        the thickening parameter for the 
                                        implicit surface.
  -P [ --project ] arg (=Newton)        defines the projection: either No or 
                                        Newton.
  -e [ --epsilon ] arg (=9.9999999999999995e-07)
                                        the maximum precision relative to the 
                                        implicit surface in the Newton 
                                        approximation of F=0.
  -n [ --max_iter ] arg (=500)          the maximum number of iteration in the 
                                        Newton approximation of F=0.
  -v [ --view ] arg (=Normal)           specifies if the surface is viewed as 
                                        is (Normal) or if places close to 
                                        singularities are highlighted 
                                        (Singular), or if unsure places should 
                                        not be displayed (Hide). @endcode

 @endcode
 
 @b Example: 
 @code
 3dImplicitSurfaceExtractorByThickening -p "x^2-y*z^2" -g 0.1 -a -2 -A 2 -v Singular
 @endcode


 You should obtain such a result:

 @image html res3dImplicitSurfaceExtractorByThickening.png "resulting visualisation."


 You could also use other implicit surfaces:
 - whitney  : x^2-y*z^2
 - 4lines   : x*y*(y-x)*(y-z*x)
 - cone     : z^2-x^2-y^2
 - simonU   : x^2-z*y^2+x^4+y^4
 - cayley3  : 4*(x^2 + y^2 + z^2) + 16*x*y*z - 1
 - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3

Some other examples (more difficult):
@code 
3dImplicitSurfaceExtractorByThickening -a -2 -A 2 -p "((y^2+z^2-1)^2-(x^2+y^2-1)^3)*(y*(x-1)^2-z*(x+1))^2" -g 0.025 -e 1e-6 -n 50000 -v Singular -t 0.5 -P Newton
3dImplicitSurfaceExtractorByThickening -a -2 -A 2 -p "(x^5-4*z^3*y^2)*((x+y)^2-(z-x)^3)" -g 0.025 -e 1e-6 -n 50000 -v Singular -t 0.05 -P Newton
 
 @endcode
 @see
 @ref 3dImplicitSurfaceExtractorByThickening.cpp

 */


template <typename CellOutputIterator, typename DigitalSurface>
void
analyzeSurface( CellOutputIterator itSure, CellOutputIterator itUnsure, DigitalSurface surface )
{
  typedef typename DigitalSurface::KSpace        KSpace;
  typedef typename DigitalSurface::Surfel        Surfel;
  typedef typename DigitalSurface::ConstIterator ConstIterator;
  typedef typename DigitalSurface::DigitalSurfaceContainer Container;
  typedef typename DigitalSurface::DigitalSurfaceTracker   Tracker;
  typedef typename KSpace::Integer               Integer;
  typedef typename KSpace::Cell                  Cell;
  const Dimension t  = KSpace::dimension - 1;
  const Container& C = surface.container();
  const KSpace& K    = surface.container().space();
  Surfel s2;
  for ( ConstIterator it = surface.begin(), itE = surface.end(); it != itE; ++it )
    {
      Surfel s   = *it;
      Cell   is  = K.unsigns( s );
      Integer s_xt = K.sKCoord( s, t );
      if ( s_xt == 0 )         *itUnsure++ = is;
      else if ( s_xt == -1 ) 
        {
          CountedPtr<Tracker> tracker( C.newTracker( s ) );
          if ( tracker->adjacent( s2, t, true ) != 0 )
            {
              Integer s2_xt = K.sKCoord( s2, t );
              Cell ic = K.uIncident( is, t, true ); 
              if ( s2_xt > 0 ) *itSure++   = ic;
              else             *itUnsure++ = ic;
            }
        }
      else if ( s_xt == 1 )
        {
          CountedPtr<Tracker> tracker( C.newTracker( s ) );
          if ( tracker->adjacent( s2, t, false ) != 0 )
            {
              Integer s2_xt = K.sKCoord( s2, t );
              Cell ic = K.uIncident( is, t, false ); 
              if ( s2_xt < 0 ) *itSure++   = ic;
              else             *itUnsure++ = ic;
            }
        }
      else cout << " " << s_xt;
    }
}

/** TODO
 * Project the point [xp] onto some zero-level of [is] with precision [eps] (meaning '|is( xp )| < eps' at the end of the process.
 * @param xp (returns) a point close to f=0.
 * @param dim the dimension of the domain.
 * @param is the implicit function.
 * @param x0 initial point.
 * @param eps the expected precision.
 */
template <typename ImplicitSurface, typename RealPoint>
RealPoint projectNewton( const ImplicitSurface & is, 
                         RealPoint x, 
                         typename RealPoint::Coordinate epsilon,
                         unsigned int max_iter )
{
  typedef typename RealPoint::Coordinate Scalar;
  RealPoint gx;
  Scalar f, g2;
  Scalar eps2 = epsilon * epsilon;
  while ( max_iter-- != 0 )
    {
      f = is( x );
      if ( abs( f ) < epsilon ) return x;
      gx = is.gradient( x );
      g2 = gx.dot( gx );
      if ( g2 < eps2 ) return x;
      x -= (f/g2) * gx;
    }
  return x;
}


template <typename CubicalComplex3, typename ImplicitShape3, 
          typename ImplicitDigitalShape3>
void projectComplex( std::vector< typename ImplicitShape3::RealPoint >& points,
                     const CubicalComplex3& complex3,
                     const ImplicitShape3& shape,
                     const ImplicitDigitalShape3& dshape,
                     double epsilon,
                     unsigned int max_iter, 
                     double max_distance )
{
  typedef typename CubicalComplex3::Cell     Cell3;
  typedef typename CubicalComplex3::Point    Point3;
  typedef typename CubicalComplex3::CellMapConstIterator CellMapConstIterator;
  typedef typename ImplicitShape3::RealPoint RealPoint3;
  typedef typename ImplicitShape3::Ring      Ring;
  points.clear();
  for ( CellMapConstIterator it = complex3.begin( 0 ), itE = complex3.end( 0 ); it != itE; ++it )
    {
      Cell3 cell    = it->first;
      Point3 dp     = complex3.space().uKCoords( cell ) - Point3::diagonal( 1 );
      RealPoint3 p  = dshape->embed( dp ) / 2.0;
      RealPoint3 q  = projectNewton( shape, p, epsilon, max_iter );
      double     d  = (p-q).norm();
      if ( d > max_distance ) q = p + (q-p)*( max_distance / d );
      points.push_back( q );
    }
}

template <typename CubicalComplex3, typename ImplicitShape3, 
          typename ImplicitDigitalShape3>
typename ImplicitShape3::Ring
getValue( const CubicalComplex3& complex3,
          const typename CubicalComplex3::Cell& cell,
          const ImplicitShape3& shape,
          const ImplicitDigitalShape3& dshape )
{
  typedef typename CubicalComplex3::Cell     Cell3;
  typedef typename CubicalComplex3::Cells    Cells3;
  typedef typename CubicalComplex3::Point    Point3;
  typedef typename ImplicitShape3::RealPoint RealPoint3;
  typedef typename ImplicitShape3::Ring      Ring;

  Point3 dp    = complex3.space().uKCoords( cell ) - Point3::diagonal( 1 );
  RealPoint3 p = dshape->embed( dp ) / 2.0;
  Ring v       = shape( p );
  return v;
}


template <typename CubicalComplex3, typename ImplicitShape3, 
          typename ImplicitDigitalShape3>
void doNotProjectComplex( std::vector< typename ImplicitShape3::RealPoint >& points,
                          const CubicalComplex3& complex3,
                          const ImplicitShape3& shape,
                          const ImplicitDigitalShape3& dshape )
{
  typedef typename CubicalComplex3::Cell     Cell3;
  typedef typename CubicalComplex3::Point    Point3;
  typedef typename CubicalComplex3::CellMapConstIterator CellMapConstIterator;
  typedef typename ImplicitShape3::RealPoint RealPoint3;
  typedef typename ImplicitShape3::Ring      Ring;
  points.clear();
  for ( CellMapConstIterator it = complex3.begin( 0 ), itE = complex3.end( 0 ); it != itE; ++it )
    {
      Cell3 cell    = it->first;
      Point3 dp     = complex3.space().uKCoords( cell );
      RealPoint3 p  = dshape->embed( dp ) / 2.0;
      points.push_back( p );
    }
}

    

int main( int argc, char** argv )
{
  typedef int                               Integer;
  typedef SpaceND<3,Integer>                Space3;
  typedef KhalimskySpaceND<3,Integer>       KSpace3;
  typedef KSpace3::Cell                     Cell3;
  typedef std::map<Cell3, CubicalCellData>  Map3;
  // typedef boost::unordered_map<Cell3, CubicalCellData>  Map3;
  typedef CubicalComplex< KSpace3, Map3 >   CC3;
  typedef Space3::Point                     Point3;
  typedef Space3::RealPoint                 RealPoint3;
  typedef Space3::RealVector                RealVector3;
  typedef RealPoint3::Coordinate            Ring;
  typedef Ring                              Scalar;
  typedef MPolynomial<3, Ring>              Polynomial3;
  typedef MPolynomialReader<3, Ring>        Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Space3>  ImplicitShape3;
  typedef GaussDigitizer< Space3, ImplicitShape3 > 
                                            ImplicitDigitalShape3;
  typedef ImplicitDigitalShape3::Domain     Domain3;
  typedef CC3::CellMapIterator              CellMapIterator;
  typedef CC3::CellMapConstIterator         CellMapConstIterator;
  typedef CC3::Cells                        Cells3;

  //-------------- parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
    ("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
    ("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " )
    ("thickness,t", po::value< double >()->default_value( 1e-2 ), "the thickening parameter for the implicit surface." )
    ("project,P", po::value< std::string >()->default_value( "Newton" ), "defines the projection: either No or Newton." )
    ("epsilon,e", po::value< double >()->default_value( 1e-6 ), "the maximum precision relative to the implicit surface in the Newton approximation of F=0." )
    ("max_iter,n", po::value< unsigned int >()->default_value( 500 ), "the maximum number of iteration in the Newton approximation of F=0." )
    ("view,v", po::value< std::string >()->default_value( "Normal" ), "specifies if the surface is viewed as is (Normal) or if places close to singularities are highlighted (Singular), or if unsure places should not be displayed (Hide)." )
    ;

  bool parseOK=true;
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  } catch ( const exception& ex ) {
    parseOK=false;
    cerr << "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if( ! parseOK || vm.count("help") || ! vm.count( "polynomial" ) )
    {
      if ( ! vm.count( "polynomial" ) ) 
        cerr << "Need parameter --polynomial / -p" << endl;

      cerr << "Usage: " << argv[0] << " -p <polynomial> [options]\n"
		<< "Computes the zero level set of the given polynomial." 
                << endl
		<< general_opt << "\n";
      cerr << "Example:\n" << endl
           << argv[0] << " -p \"x^2-y*z^2\" -g 0.1 -a -2 -A 2 -v Singular" << endl
           << " - whitney  : x^2-y*z^2" << endl
           << " - 4lines   : x*y*(y-x)*(y-z*x)" << endl
           << " - cone     : z^2-x^2-y^2" << endl
           << " - simonU   : x^2-z*y^2+x^4+y^4" << endl
           << " - cayley3  : 4*(x^2 + y^2 + z^2) + 16*x*y*z - 1" << endl
           << " - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3" << endl << endl;
      cerr << "Some other examples (more difficult):" << endl
           << argv[0] << " -a -2 -A 2 -p \"((y^2+z^2-1)^2-(x^2+y^2-1)^3)*(y*(x-1)^2-z*(x+1))^2\" -g 0.025 -e 1e-6 -n 50000 -v Singular -t 0.5 -P Newton" << endl
           << argv[0] << " -a -2 -A 2 -p \"(x^5-4*z^3*y^2)*((x+y)^2-(z-x)^3)\" -g 0.025 -e 1e-6 -n 50000 -v Singular -t 0.05 -P Newton" << endl;
      return 0;
    }

  //-------------- read polynomial and creating 3d implicit fct -----------------
  trace.beginBlock( "Reading polynomial and creating 3D implicit function" );
  string poly_str = vm[ "polynomial" ].as<string>();
  Polynomial3 poly;
  Polynomial3Reader reader;
  string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
    {
      trace.error() << "ERROR reading polynomial: I read only <" << poly_str.substr( 0, iter - poly_str.begin() )
                    << ">, and I built P=" << poly << std::endl;
      return 2;
    }
  // Creating implicit shape and storing it with smart pointer for automatic deallocation.
  CountedPtr<ImplicitShape3> shape( new ImplicitShape3( poly ) );

  Ring min_x = vm[ "minAABB" ].as<double>();
  Ring max_x = vm[ "maxAABB" ].as<double>();
  Ring h     = vm[ "gridstep" ].as<double>();
  RealPoint3 p1( min_x, min_x, min_x );
  RealPoint3 p2( max_x, max_x, max_x );
  // Creating digitized shape and storing it with smart pointer for automatic deallocation.
  CountedPtr<ImplicitDigitalShape3> dshape( new ImplicitDigitalShape3() );
  dshape->attach( *shape );
  dshape->init( p1, p2, RealVector3( h, h, h ) );
  Domain3 domain3 = dshape->getDomain();
  KSpace3 K3;
  K3.init( domain3.lowerBound(), domain3.upperBound(), true );
  trace.info() << "- domain is " << domain3 << std::endl;
  trace.endBlock();

  //-------------- read polynomial and creating 3d implicit fct -----------------
  trace.beginBlock( "Extracting thickened isosurface [-t,+t] of 3D polynomial. " );
  CubicalCellData unsure_data( 0 );
  CubicalCellData sure_data( CC3::FIXED );
  Ring t     = vm[ "thickness" ].as<double>();
  CC3  complex3( K3 );
  for ( Domain3::ConstIterator it = domain3.begin(), itE = domain3.end(); it != itE; ++it )
    {
      Cell3 spel    = K3.uSpel( *it );
      // RealPoint3 px = dshape->embed( *it );
      RealPoint3 px = dshape->embed( K3.uKCoords( spel ) - Point3::diagonal( 1 ) ) / 2.0;
      Ring s        = (*shape)( px );
      if ( (-t <= s ) && ( s <= t ) ) complex3.insertCell( spel, unsure_data );
    }
  trace.info() << "-    K[-t,+t]        = " << complex3 << endl;
  complex3.close();
  trace.info() << "- Cl K[-t,+t]        = " << complex3 << endl;
  std::vector<Cell3> separating_cells;
  std::back_insert_iterator< std::vector<Cell3> > outItSurface( separating_cells );
  Surfaces<KSpace3>::uWriteBoundary( outItSurface,
                                     K3, *dshape, domain3.lowerBound(), domain3.upperBound() );
  trace.info() << "- separating S       = " << separating_cells.size() << " 2-cells." << endl;
  complex3.insertCells( separating_cells.begin(), separating_cells.end(), sure_data );
  for ( std::vector<Cell3>::const_iterator it = separating_cells.begin(), itE = separating_cells.end(); it != itE; ++it )
    {
      Cells3 bdry = K3.uFaces( *it );
      for ( Cells3::const_iterator itBdry = bdry.begin(), itBdryE = bdry.end(); itBdry != itBdryE; ++itBdry )
        complex3.insertCell( *itBdry, sure_data );
    }
  separating_cells.clear();
  trace.info() << "- Cl K[-t,+t] + Cl S = " << complex3 << endl;
  trace.endBlock();

  //-------------- Get boundary and inner cells  --------------------------------
  trace.beginBlock( "Get boundary and inner cells. " );
  std::vector<Cell3> inner;
  std::vector<Cell3> bdry;
  functions::filterCellsWithinBounds
    ( complex3, K3.uKCoords( K3.lowerCell() ), K3.uKCoords( K3.upperCell() ),
       std::back_inserter( bdry ), std::back_inserter( inner ) );
  trace.info() << "- there are " << inner.size() << " inner cells." << endl;
  trace.info() << "- there are " << bdry.size() << " boundary cells." << endl;
  trace.endBlock();

  //-------------- Compute priority function  -----------------------------------
  trace.beginBlock( "Compute priority function. " );
  Dimension d = complex3.dim();
  for ( Dimension i = 0; i <= d; ++i )
    {
      for ( CellMapIterator it = complex3.begin( i ), itE = complex3.end( i ); it != itE; ++it )
        {
          Cell3 cell   = it->first;
          Ring v       = getValue( complex3, cell, *shape, dshape );
          v = abs( 10000.0*v );
          if ( v > 10000000.0 ) v = 10000000.0;
          it->second.data &= ~CC3::VALUE;
          it->second.data |= (DGtal::uint32_t) floor( v );
          // std::cout << " " << it->second.data;
        }
    }
  trace.endBlock();

  //-------------- Collapse boundary -------------------------------------------
  trace.beginBlock( "Collapse boundary. " );
  typename CC3::DefaultCellMapIteratorPriority priority;
  CC3 bdry_complex3( K3 );
  for ( std::vector<Cell3>::const_iterator it = bdry.begin(), itE = bdry.end(); it != itE; ++it )
    {
      Cell3 cell = *it;
      Dimension d = K3.uDim( cell );
      CellMapConstIterator cmIt = complex3.findCell( d, cell );
      bdry_complex3.insertCell( d, cell, cmIt->second );
    }
  trace.info() << "- [before collapse] K_bdry =" << bdry_complex3 << endl;
  functions::collapse( bdry_complex3, bdry.begin(), bdry.end(), priority, true, true, false );
  trace.info() << "- [after collapse]  K_bdry =" << bdry_complex3 << endl;
  for ( std::vector<Cell3>::const_iterator it = bdry.begin(), itE = bdry.end(); it != itE; ++it )
    {
      Cell3 cell  = *it;
      Dimension d = K3.uDim( cell );
      CellMapConstIterator cmIt = bdry_complex3.findCell( d, cell );
      if ( cmIt != bdry_complex3.end( d ) ) {
        CellMapIterator cmIt2 = complex3.findCell( d, cell );
        cmIt2->second = sure_data;
      }
    }
  trace.endBlock();

  //-------------- Collapse all -------------------------------------------
  trace.beginBlock( "Collapse all. " );
  std::copy( bdry.begin(), bdry.end(), std::back_inserter( inner ) );
  functions::collapse( complex3, inner.begin(), inner.end(), priority, true, true, true );
  trace.info() << "- K = " << complex3 << endl;
  trace.endBlock();

  //-------------- Project complex onto surface --------------------------------
  trace.beginBlock( "Project complex onto surface. " );
  std::string project   = vm[ "project" ].as<std::string>();
  double epsilon        = vm[ "epsilon" ].as<double>();
  unsigned int max_iter = vm[ "max_iter" ].as<unsigned int>();
  std::vector<RealPoint3> points;
  if ( project == "Newton" )
    projectComplex( points, complex3, *shape, dshape, epsilon, max_iter, h * sqrt(3.0));
  else
    doNotProjectComplex( points, complex3, *shape, dshape );
  trace.endBlock();

  //-------------- Create Mesh -------------------------------------------
  trace.beginBlock( "Create Mesh. " );
  std::string view = vm[ "view" ].as<std::string>();
  bool highlight = ( view == "Singular" );
  bool hide      = ( view == "Hide" );
  Mesh<RealPoint3> mesh( true );
  std::map<Cell3,unsigned int> indices;
  int idx = 0;
  for ( CellMapConstIterator it = complex3.begin( 0 ), itE = complex3.end( 0 ); it != itE; ++it, ++idx )
    {
      Cell3 cell = it->first;
      indices[ cell ] = idx;
      mesh.addVertex( points[ idx ] );
    }
  for ( CellMapConstIterator it = complex3.begin( 2 ), itE = complex3.end( 2 ); it != itE; ++it )
    {
      Cell3 cell = it->first;
      bool fixed = it->second.data & CC3::FIXED;
      Cells3 bdry = complex3.cellBoundary( cell, true );
      std::vector<unsigned int> face_idx;
      for ( Cells3::const_iterator itC = bdry.begin(), itCE = bdry.end(); itC != itCE; ++itC )
        {
          if ( complex3.dim( *itC ) == 0 )
            face_idx.push_back( indices[ *itC ] );
        }
      if ( ( ! fixed ) && hide ) continue;
      Color color = highlight
        ? ( fixed ? Color::White : Color(128,255,128) )
        : Color::White;
      RealVector3 diag03 = points[ face_idx[ 0 ] ] - points[ face_idx[ 3 ] ];
      RealVector3 diag12 = points[ face_idx[ 1 ] ] - points[ face_idx[ 2 ] ];
      if ( diag03.dot( diag03 ) <= diag12.dot( diag12 ) )
        {
          mesh.addTriangularFace( face_idx[ 0 ], face_idx[ 1 ], face_idx[ 3 ], color );
          mesh.addTriangularFace( face_idx[ 0 ], face_idx[ 3 ], face_idx[ 2 ], color );
        }
      else
        {
          mesh.addTriangularFace( face_idx[ 0 ], face_idx[ 1 ], face_idx[ 2 ], color );
          mesh.addTriangularFace( face_idx[ 1 ], face_idx[ 3 ], face_idx[ 2 ], color );
        }
      //mesh.addQuadFace( face_idx[ 0 ], face_idx[ 1 ], face_idx[ 3 ], face_idx[ 2 ], color );
     }
  trace.endBlock();

  //-------------- View surface -------------------------------------------
  QApplication application(argc,argv);
  Viewer3D<Space3,KSpace3> viewer( K3 );
  viewer.setWindowTitle("Implicit surface viewer by thickening");
  viewer.show();
  viewer << mesh;
  // Display lines that are not in the mesh.
  for ( CellMapConstIterator it = complex3.begin( 1 ), itE = complex3.end( 1 ); it != itE; ++it )
    {
      Cell3 cell  = it->first;
      bool fixed  = it->second.data & CC3::FIXED;
      std::vector<Cell3> dummy;
      std::back_insert_iterator< std::vector<Cell3> > outIt( dummy );
      complex3.directCoFaces( outIt, cell );
      if ( ! dummy.empty() )     continue;

      Cells3 bdry = complex3.cellBoundary( cell, true );
      Cell3 v0    = *(bdry.begin() );
      Cell3 v1    = *(bdry.begin() + 1);
      if ( ( ! fixed ) && hide ) continue;
      Color color = highlight
        ? ( fixed ? Color::White : Color(128,255,128) )
        : Color::White;
      viewer.setLineColor( color );
      viewer.addLine( points[ indices[ v0 ] ], points[ indices[ v1 ] ], h/2.0 );
    }
  viewer << Viewer3D<Space3,KSpace3>::updateDisplay;
  return application.exec();
}
///////////////////////////////////////////////////////////////////////////////
