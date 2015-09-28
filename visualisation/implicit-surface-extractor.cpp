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
 * @file cubical-complex-collapse.cpp
 * @ingroup Examples
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/08/28
 *
 * An example file named cubical-complex-collapse.cpp.
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
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/shapes/GaussDigitizer.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;


/**
 * This 4D extension of a 3D implicit function transforms f(x,y,z) as:
 * \f$ F(x,y,z,t) = f(x,y,z) - t | \nabla f(x,y,z) | \f$.
 */
template <typename TSpace3, typename TSpace4, typename TPolynomial3>
struct ImplicitSurface4DExtension {
  typedef TSpace3                          Space3;
  typedef TSpace4                          Space4;
  typedef TPolynomial3                     Polynomial3;
  typedef typename Space3::RealPoint       RealPoint3;
  typedef typename Space3::RealVector      RealVector3;
  typedef typename Space3::Integer         Integer;
  typedef typename RealPoint3::Coordinate  Ring;
  typedef typename Space4::RealPoint       RealPoint4;
  typedef typename Space4::RealVector      RealVector4;
  typedef RealPoint4                       RealPoint;
  typedef RealVector4                      RealVector;

  ImplicitSurface4DExtension( Clone<Polynomial3> poly )
    : f( poly )
  {
    fx  = derivative<0>( f );
    fy  = derivative<1>( f );
    fz  = derivative<2>( f );
    fxx = derivative<0>( fx );
    fxy = derivative<0>( fy );
    fxz = derivative<0>( fz );
    fyy = derivative<1>( fy );
    fyz = derivative<1>( fz );
    fzz = derivative<2>( fz );
  }

  /**
     @param aPoint any point in the Euclidean space.
     @return the value of the polynomial at \a aPoint.
  */
  Ring operator()(const RealPoint &aPoint) const
  {
    Ring x = aPoint[ 0 ];
    Ring y = aPoint[ 1 ];
    Ring z = aPoint[ 2 ];
    Ring t = aPoint[ 3 ];
    RealVector3 grad( fx(x)(y)(z), fy(x)(y)(z), fz(x)(y)(z) );
    return f(x)(y)(z) - t * grad.norm();
  }

  /**
     @param aPoint any point in the Euclidean space.
     @return 'true' if the polynomial value is < 0.
  */
  bool isInside(const RealPoint &aPoint) const
  {
    return this->operator()( aPoint ) <= NumberTraits<Ring>::ZERO;
  }

  /**
     @param aPoint any point in the Euclidean space.
     
     @return INSIDE if the polynomial value is < 0, OUTSIDE if > 0,
     ON otherwise.
  */
  Orientation orientation(const RealPoint &aPoint) const
  {
    Ring v = this->operator()( aPoint );
    if ( v > NumberTraits<Ring>::ZERO )      return OUTSIDE;
    else if ( v < NumberTraits<Ring>::ZERO ) return INSIDE;
    else                                     return ON;
  }

  /**
     @param aPoint any point in the Euclidean space.
     @return the gradient vector of the polynomial at \a aPoint.
  */
  inline
  RealVector gradient( const RealPoint &aPoint ) const
  {
    Ring x = aPoint[ 0 ];
    Ring y = aPoint[ 1 ];
    Ring z = aPoint[ 2 ];
    Ring t = aPoint[ 3 ];
    RealVector3 grad( fx(x)(y)(z), fy(x)(y)(z), fz(x)(y)(z) );
    Ring ngrad = grad.norm();
    grad /= ngrad;
    RealVector d_gradx( fxx(x)(y)(z), fxy(x)(y)(z),fxz(x)(y)(z) );
    RealVector d_grady( fxy(x)(y)(z), fyy(x)(y)(z),fyz(x)(y)(z) );
    RealVector d_gradz( fxz(x)(y)(z), fyz(x)(y)(z),fzz(x)(y)(z) );
    return RealVector( grad[ 0 ] - t * d_gradx.dot( grad ),
                       grad[ 1 ] - t * d_grady.dot( grad ),
                       grad[ 2 ] - t * d_gradz.dot( grad ),
                       -ngrad );
                      
  }

private:
  Polynomial3 f;
  Polynomial3 fx;
  Polynomial3 fy;
  Polynomial3 fz;
  Polynomial3 fxx;
  Polynomial3 fxy;
  Polynomial3 fxz;
  Polynomial3 fyy;
  Polynomial3 fyz;
  Polynomial3 fzz;
};

/**
 * Computes the cells of the given complex that lies on the
 * boundary or inside the given bounds.
 */
template <typename CubicalComplex, 
          typename CellOutputIterator>
void 
getCellsWithinBounds( CellOutputIterator itBdry, CellOutputIterator itInner,
                   const CubicalComplex& K, 
                   const typename CubicalComplex::Point& kLow,  
                   const typename CubicalComplex::Point& kUp )
{
  typedef typename CubicalComplex::Cell                 Cell;
  typedef typename CubicalComplex::Point                Point;
  typedef typename CubicalComplex::CellMapConstIterator CellMapConstIterator;
  Dimension d = K.dim();
  for ( Dimension i = 0; i <= d; ++i )
    {
      for ( CellMapConstIterator it = K.begin( i ), itE = K.end( i ); it != itE; ++it )
        {
          Cell cell = it->first;
          Point kCell = K.space().uKCoords( cell );
          if ( ( kLow.inf( kCell ) == kLow ) && ( kUp.sup( kCell ) == kUp ) )
            { // Inside or on boundary.
              bool bdry = false;
              for ( Dimension j = 0; j < Point::dimension - 1; ++j )
                {
                  if ( ( kCell[ j ] == kLow[ j ] ) || ( kCell[ j ] == kUp[ j ] ) )
                    {
                      bdry = true;
                      break;
                    }
                }
              if ( bdry ) *itBdry++  = cell;
              else        *itInner++ = cell;
            }
        }
    }
}

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

int main( int argc, char** argv )
{
  typedef int                               Integer;
  typedef SpaceND<3,Integer>                Space3;
  typedef KhalimskySpaceND<3,Integer>       KSpace3;
  typedef KSpace3::Cell                     Cell3;
  typedef std::map<Cell3, CubicalCellData>  Map3;
  typedef CubicalComplex< KSpace3, Map3 >   CC3;
  typedef Space3::Point                     Point3;
  typedef Space3::RealPoint                 RealPoint3;
  typedef RealPoint3::Coordinate            Ring;
  typedef Ring                              Scalar;
  typedef MPolynomial<3, Ring>              Polynomial3;
  typedef MPolynomialReader<3, Ring>        Polynomial3Reader;
  typedef SpaceND<4,Integer>                Space4;
  typedef KhalimskySpaceND<4,Integer>       KSpace4;
  typedef KSpace4::Cell                     Cell4;
  typedef KSpace4::SCell                    SCell4;
  typedef Space4::Point                     Point4;
  typedef Space4::RealPoint                 RealPoint4;
  typedef ImplicitSurface4DExtension<Space3,Space4,Polynomial3>
                                            ImplicitShape4;
  typedef GaussDigitizer< Space4, ImplicitShape4 > 
                                            ImplicitDigitalShape4;
  typedef ImplicitDigitalShape4::Domain     Domain4;
  typedef std::map<Cell4, CubicalCellData>  Map4;
  typedef CubicalComplex< KSpace4, Map4 >   CC4;
  typedef CC4::CellMapIterator              CellMapIterator;
  typedef CC4::CellMapConstIterator         CellMapConstIterator;

  //-------------- parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
    ("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
    ("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " )
    ("timestep,t", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization in the 4th dimension. " )
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
           << "./implicit-surface-extractor -p \"-0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3\"" << endl
           << " - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3" << endl << endl;
      return 0;
    }

  //-------------- read polynomial and creating 4d implicit fct -----------------
  trace.beginBlock( "Reading polynomial and creating 4D implicit function" );
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
  CountedPtr<ImplicitShape4> shape( new ImplicitShape4( poly ) ); 

  Ring min_x = vm[ "minAABB" ].as<double>();
  Ring max_x = vm[ "maxAABB" ].as<double>();
  Ring h     = vm[ "gridstep" ].as<double>();
  Ring t     = vm[ "timestep" ].as<double>();
  RealPoint4 p1( min_x, min_x, min_x, -t*h );
  RealPoint4 p2( max_x, max_x, max_x,  0 );
  // Creating digitized shape and storing it with smart pointer for automatic deallocation.
  CountedPtr<ImplicitDigitalShape4> dshape( new ImplicitDigitalShape4() );
  dshape->attach( *shape );
  dshape->init( p1, p2, h );
  Domain4 domain4 = dshape->getDomain();
  KSpace4 K4;
  K4.init( domain4.lowerBound(), domain4.upperBound(), true );
  trace.info() << "- domain is " << domain4 << std::endl;
  trace.endBlock();

  //-------------- read polynomial and creating 4d implicit fct -----------------
  trace.beginBlock( "Extracting isosurface 0 of 4D polynomial. " );
  CC4 complex4( K4 );
  std::vector<SCell4> surface_cells;
  std::back_insert_iterator< std::vector<SCell4> > outItSurface( surface_cells );
  Surfaces<KSpace4>::sWriteBoundary( outItSurface,
                                     K4, *dshape, domain4.lowerBound(), domain4.upperBound() );
  trace.info() << "- 4D surface has " << surface_cells.size() << " 3-cells." << endl;
  typedef SetOfSurfels<KSpace4> SurfaceContainer;
  SurfaceContainer* ptrF0 = new SurfaceContainer( K4, SurfelAdjacency< KSpace4::dimension >( false ) );
  ptrF0->surfelSet().insert( surface_cells.begin(), surface_cells.end() );
  DigitalSurface<SurfaceContainer> surface( ptrF0 );
  std::vector<Cell4> sure_cells;
  std::vector<Cell4> unsure_cells;
  analyzeSurface( std::back_inserter( sure_cells ), std::back_inserter( unsure_cells ),
                  surface );
  trace.info() << "- sure cells   = " << sure_cells.size() << endl;
  trace.info() << "- unsure cells = " << unsure_cells.size() << endl;
  CubicalCellData unsure_data( 0 );
  CubicalCellData sure_data( CC4::FIXED );
  complex4.insertCells( sure_cells.begin(), sure_cells.end(), sure_data );
  complex4.insertCells( unsure_cells.begin(), unsure_cells.end(), unsure_data );
  surface_cells.clear();
  trace.endBlock();

  //-------------- Get boundary and inner cells  --------------------------------
  trace.beginBlock( "Get boundary and inner cells. " );
  std::vector<Cell4> inner;
  std::vector<Cell4> bdry;
  complex4.close();
  trace.info() << "- K=" << complex4 << endl;
  getCellsWithinBounds( std::back_inserter( bdry ), std::back_inserter( inner ),
                        complex4, K4.uKCoords( K4.lowerCell() ), K4.uKCoords( K4.upperCell() ) );
  trace.info() << "- there are " << inner.size() << " inner cells." << endl;
  trace.info() << "- there are " << bdry.size() << " boundary cells." << endl;
  trace.endBlock();

  //-------------- Collapse boundary -------------------------------------------
  trace.beginBlock( "Collapse boundary. " );
  typename CC4::DefaultCellMapIteratorPriority priority;
  CC4 bdry_complex4( K4 );
  bdry_complex4.insertCells( bdry.begin(), bdry.end() );
  trace.info() << "- [before collapse] K_bdry =" << bdry_complex4 << endl;
  bdry_complex4.collapse( bdry.begin(), bdry.end(), priority, true, true, true );
  trace.info() << "- [after collapse]  K_bdry =" << bdry_complex4 << endl;
  for ( std::vector<Cell4>::const_iterator it = bdry.begin(), itE = bdry.end(); it != itE; ++it )
    {
      Cell4 cell = *it;
      Dimension d = K4.uDim( cell );
      CellMapConstIterator cmIt = bdry_complex4.find( d, cell );
      if ( cmIt != bdry_complex4.end( d ) ) {
        CellMapIterator cmIt2 = complex4.find( d, cell );
        cmIt2->second = sure_data;
      }
    }
  trace.endBlock();

  //-------------- Collapse all -------------------------------------------
  trace.beginBlock( "Collapse all. " );
  std::copy( bdry.begin(), bdry.end(), std::back_inserter( inner ) );
  //typename CC4::DefaultCellMapIteratorPriority priority;
  complex4.collapse( inner.begin(), inner.end(), priority, false, false, true );
  trace.info() << "- K=" << complex4 << endl;
  trace.endBlock();

  QApplication application(argc,argv);
  Point4 low4 = K4.lowerBound();
  Point4 up4  = K4.upperBound();
  KSpace3 K3;
  K3.init( Point3( low4[ 0 ], low4[ 1 ], low4[ 2 ] ),
           Point3( up4 [ 0 ], up4 [ 1 ], up4 [ 2 ] ), true );
  Viewer3D<Space3,KSpace3> viewer( K3 );
  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();
  Dimension d = complex4.dim();
  for ( Dimension i = 0; i <= d; ++i )
    {
      for ( CellMapConstIterator it = complex4.begin( i ), itE = complex4.end( i ); it != itE; ++it )
        {
          bool fixed = it->second.data & CC4::FIXED;
          Cell4 cell = it->first;
          Point4 kcell = K4.uKCoords( cell ); 
          Cell3 proj_cell = K3.uCell( Point3( kcell[ 0 ], kcell[ 1 ], kcell[ 2 ] ) );
          viewer.setFillColor( fixed ? Color::Red : Color::White );
          viewer << proj_cell;
        }
    }
  viewer << Viewer3D<Space3,KSpace3>::updateDisplay;
  return application.exec();


  // trace.beginBlock( "Creating Cubical Complex" );
  // KSpace K;
  // K.init( Point( 0,0,0 ), Point( 512,512,512 ), true );
  // CC complex( K );
  // Integer m = 10;
  // std::vector<Cell> S;
  // for ( Integer x = 0; x <= m; ++x )
  //   for ( Integer y = 0; y <= m; ++y )
  //     for ( Integer z = 0; z <= m; ++z )
  //       {
  //         Point k1 = Point( x, y, z ); 
  //         S.push_back( K.uCell( k1 ) );
  //         double d1 = Point::diagonal( 1 ).dot( k1 ) / sqrt( (double) KSpace::dimension );
  //         RealPoint v1( k1[ 0 ] - d1 * k1[ 0 ], k1[ 1 ] - d1 * k1[ 1 ], k1[ 2 ] - d1 * k1[ 2 ] );
  //         double n1 = v1.dot( v1 );
  //         bool fixed = ( ( x == 0 ) && ( y == 0 ) && ( z == 0 ) )
  //           || ( ( x == 0 ) && ( y == m ) && ( z == 0 ) )
  //           || ( ( x == m ) && ( y == 0 ) && ( z == 0 ) )
  //           || ( ( x == m ) && ( y == m ) && ( z == 0 ) )
  //           || ( ( x == m/3 ) && ( y == 2*m/3 ) && ( z == 2*m/3 ) )
  //           || ( ( x == 0 ) && ( y == 0 ) && ( z == m ) )
  //           || ( ( x == 0 ) && ( y == m ) && ( z == m ) )
  //           || ( ( x == m ) && ( y == 0 ) && ( z == m ) )
  //           || ( ( x == m ) && ( y == m ) && ( z == m ) );
  //         complex.insertCell( S.back(), 
  //                             fixed ? CC::FIXED 
  //                             : (uint32_t) floor(64.0 * sqrt( n1 ) ) // This is the priority for collapse 
  //                             );
  //       }
  // //complex.close();
  // trace.info() << "After close: " << complex << std::endl;
  // trace.endBlock();

  // // for 3D display with Viewer3D
  // QApplication application(argc,argv);
  // typedef Viewer3D<Space, KSpace> MyViewer;

  // {
  //   MyViewer viewer(K);
  //   viewer.show();
  //   typedef CC::CellMapConstIterator CellMapConstIterator;
  //   for ( Dimension d = 0; d <= 3; ++d )
  //     for ( CellMapConstIterator it = complex.begin( d ), itE = complex.end( d );
  //           it != itE; ++it )
  //       {
  //         viewer << it->first;
  //       }
  //   viewer<< MyViewer::updateDisplay;
  //   application.exec();
  // }
  
  // trace.beginBlock( "Collapsing complex" );
  // CC::DefaultCellMapIteratorPriority P;
  // //DiagonalPriority<CC> P( complex );
  // complex.collapse( S.begin(), S.end(), P, true, true );
  // trace.info() << "After collapse: " << complex << std::endl;
  // trace.endBlock();

  // {
  //   MyViewer viewer(K);
  //   viewer.show();
  //   typedef CC::CellMapConstIterator CellMapConstIterator;
  //   for ( Dimension d = 0; d <= 3; ++d )
  //     for ( CellMapConstIterator it = complex.begin( d ), itE = complex.end( d );
  //           it != itE; ++it )
  //       {
  //         viewer << it->first;
  //       }
  //   viewer<< MyViewer::updateDisplay;
  //   return application.exec();
  // }
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
