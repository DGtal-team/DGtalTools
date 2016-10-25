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
 * @file 3dImplicitSurfaceExtractorBy4DExtension.cpp
 * @ingroup Visualisation
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/08/28
 *
 * A tool to visualize 3D implicit surface by viewing it as 4D
 * hypersurface to detect zeroes and then shrink it with cubical
 * complex collapses.
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
 @page Doc3dImplicitSurfaceExtractorBy4DExtension 3dImplicitSurfaceExtractorBy4DExtension
 
 @brief Computes the zero level set of the given polynomial.

 @b Usage:  3dImplicitSurfaceExtractorBy4DExtension [options] input

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
  -t [ --timestep ] arg (=9.9999999999999995e-07)
                                        the gridstep that defines the 
                                        digitization in the 4th dimension 
                                        (small is generally a good idea, 
                                        default is 1e-6). 
  -P [ --project ] arg (=Newton)        defines the projection: either No or 
                                        Newton.
  -e [ --epsilon ] arg (=9.9999999999999995e-08)
                                        the maximum precision relative to the 
                                        implicit surface.
  -n [ --max_iter ] arg (=500)          the maximum number of iteration in the 
                                        Newton approximation of F=0.
  -v [ --view ] arg (=Normal)           specifies if the surface is viewed as 
                                        is (Normal) or if places close to 
                                        singularities are highlighted 
                                        (Singular).
 @endcode


 @b Example: 

 @code
   3dImplicitSurfaceExtractorBy4DExtension -p "-0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3" -g 0.06125 -a -2 -A 2 -v Singular -t 0.02
 @endcode
 You could also use other implicit surfaces:
 - whitney  : x^2-y*z^2
 - 4lines   : x*y*(y-x)*(y-z*x)
 - cone     : z^2-x^2-y^2
 - simonU   : x^2-z*y^2+x^4+y^4
 - cayley3  : 4*(x^2 + y^2 + z^2) + 16*x*y*z - 1
 - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3


 You should obtain such a result:

 @image html res3dImplicitSurfaceExtractorBy4DExtension.png "resulting visualisation."
 

 @see
 @ref 3dImplicitSurfaceExtractorBy4DExtension.cpp

 */



/**
 * This 4D extension of a 3D implicit function transforms f(x,y,z) as:
 * \f$ F(x,y,z,t) = f(x,y,z) - t | \nabla f(x,y,z) | \f$.
 */
template <typename TSpace3, typename TSpace4, typename TPolynomial3>
struct ImplicitSurface4DProductExtension {
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

  ImplicitSurface4DProductExtension( Clone<Polynomial3> poly )
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
 * This 4D extension of a 3D implicit function transforms f(x,y,z) as:
 * \f$ F(x,y,z,t) = f(x,y,z) - t \f$.
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
    // RealVector3 grad( fx(x)(y)(z), fy(x)(y)(z), fz(x)(y)(z) );
    return f(x)(y)(z) - t;
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
    return RealVector( fx(x)(y)(z), fy(x)(y)(z), fz(x)(y)(z), -1.0 );
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


template <typename CellOutputIterator, typename DigitalSurface>
void
analyzeSurface( CellOutputIterator itSure, CellOutputIterator itUnsure, 
                DigitalSurface surface )
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


template <typename CubicalComplex4, typename ImplicitShape4, 
          typename ImplicitDigitalShape4, typename ImplicitShape3>
void projectComplex( std::vector< typename ImplicitShape3::RealPoint >& points,
                     const CubicalComplex4& complex4,
                     const ImplicitShape4& shape,
                     const ImplicitDigitalShape4& dshape,
                     const ImplicitShape3& shape3,
                     double epsilon,
                     unsigned int max_iter )
{
  typedef typename CubicalComplex4::Cell     Cell4;
  typedef typename CubicalComplex4::Point    Point4;
  typedef typename CubicalComplex4::CellMapConstIterator CellMapConstIterator;
  typedef typename ImplicitShape4::RealPoint RealPoint4;
  typedef typename ImplicitShape4::Ring      Ring;
  typedef typename ImplicitShape3::RealPoint RealPoint3;
  points.clear();
  for ( CellMapConstIterator it = complex4.begin( 0 ), itE = complex4.end( 0 ); it != itE; ++it )
    {
      Cell4 cell    = it->first;
      Point4 dp     = complex4.space().uCoords( cell );
      RealPoint4 p  = dshape->embed( dp );
      RealPoint3 p3 = RealPoint3( p[ 0 ], p[ 1 ], p[ 2 ] );
      RealPoint3 q  = projectNewton( shape3, p3, epsilon, max_iter );
      points.push_back( q );
    }
}

template <typename CubicalComplex4, typename ImplicitShape4, 
          typename ImplicitDigitalShape4, typename ImplicitShape3>
void doNotProjectComplex( std::vector< typename ImplicitShape3::RealPoint >& points,
                          const CubicalComplex4& complex4,
                          const ImplicitShape4& shape,
                          const ImplicitDigitalShape4& dshape,
                          const ImplicitShape3& shape3 )
{
  typedef typename CubicalComplex4::Cell     Cell4;
  typedef typename CubicalComplex4::Point    Point4;
  typedef typename CubicalComplex4::CellMapConstIterator CellMapConstIterator;
  typedef typename ImplicitShape4::RealPoint RealPoint4;
  typedef typename ImplicitShape4::Ring      Ring;
  typedef typename ImplicitShape3::RealPoint RealPoint3;
  for ( CellMapConstIterator it = complex4.begin( 0 ), itE = complex4.end( 0 ); it != itE; ++it )
    {
      Cell4 cell    = it->first;
      Point4 dp     = complex4.space().uCoords( cell );
      RealPoint4 p  = dshape->embed( dp );
      RealPoint3 p3 = RealPoint3( p[ 0 ], p[ 1 ], p[ 2 ] );
      points.push_back( p3 );
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
  typedef ImplicitPolynomial3Shape<Space3>  ImplicitShape3;
  typedef SpaceND<4,Integer>                Space4;
  typedef KhalimskySpaceND<4,Integer>       KSpace4;
  typedef KSpace4::Cell                     Cell4;
  typedef KSpace4::Cells                    Cells4;
  typedef KSpace4::SCell                    SCell4;
  typedef Space4::Point                     Point4;
  typedef Space4::RealPoint                 RealPoint4;
  typedef Space4::RealVector                RealVector4;
  typedef ImplicitSurface4DProductExtension<Space3,Space4,Polynomial3>
                                            ImplicitShape4;
  typedef GaussDigitizer< Space4, ImplicitShape4 > 
                                            ImplicitDigitalShape4;
  typedef ImplicitDigitalShape4::Domain     Domain4;
  typedef std::map<Cell4, CubicalCellData>  Map4;
  typedef CubicalComplex< KSpace4, Map4 >   CC4;
  typedef CC4::CellMapIterator              CellMapIterator;
  typedef CC4::CellMapConstIterator         CellMapConstIterator;
  typedef CC4::Cells                        Cells4;

  //-------------- parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
    ("minAABB,a",  po::value<double>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
    ("maxAABB,A",  po::value<double>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " )
    ("timestep,t", po::value< double >()->default_value( 0.000001 ), "the gridstep that defines the digitization in the 4th dimension (small is generally a good idea, default is 1e-6). " )
    ("project,P", po::value< std::string >()->default_value( "Newton" ), "defines the projection: either No or Newton." )
    ("epsilon,e", po::value< double >()->default_value( 0.0000001 ), "the maximum precision relative to the implicit surface." )
    ("max_iter,n", po::value< unsigned int >()->default_value( 500 ), "the maximum number of iteration in the Newton approximation of F=0." )
    ("view,v", po::value< std::string >()->default_value( "Normal" ), "specifies if the surface is viewed as is (Normal) or if places close to singularities are highlighted (Singular)." )
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
           << argv[ 0 ] << " -p \"-0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3\" -g 0.06125 -a -2 -A 2 -v Singular -t 0.02" << endl
           << " - whitney  : x^2-y*z^2" << endl
           << " - 4lines   : x*y*(y-x)*(y-z*x)" << endl
           << " - cone     : z^2-x^2-y^2" << endl
           << " - simonU   : x^2-z*y^2+x^4+y^4" << endl
           << " - cayley3  : 4*(x^2 + y^2 + z^2) + 16*x*y*z - 1" << endl
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
  ImplicitShape3 shape3( poly );
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
  dshape->init( p1, p2, RealVector4( h, h, h, t ) );
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
  functions::filterCellsWithinBounds
    ( complex4, K4.uKCoords( K4.lowerCell() ), K4.uKCoords( K4.upperCell() ),
      std::back_inserter( bdry ), std::back_inserter( inner ) );
  trace.info() << "- there are " << inner.size() << " inner cells." << endl;
  trace.info() << "- there are " << bdry.size() << " boundary cells." << endl;
  trace.endBlock();

  //-------------- Compute priority function  -----------------------------------
  trace.beginBlock( "Compute priority function. " );
  Dimension d = complex4.dim();
  for ( Dimension i = 0; i <= d; ++i )
    {
      for ( CellMapIterator it = complex4.begin( i ), itE = complex4.end( i ); it != itE; ++it )
        {
          Cell4 cell   = it->first;
          Point4 dp    = K4.uCoords( cell );
          RealPoint4 p = dshape->embed( dp );
          Ring v       = (*shape)( p );
          v = abs( 1000.0*v );
          if ( v > 1000000.0 ) v = 1000000.0;
          it->second.data &= ~CC4::VALUE;
          it->second.data |= (DGtal::uint32_t) floor( v );
          // std::cout << " " << it->second.data;
        }
    }
  trace.endBlock();

  //-------------- Collapse boundary -------------------------------------------
  trace.beginBlock( "Collapse boundary. " );
  typename CC4::DefaultCellMapIteratorPriority priority;
  CC4 bdry_complex4( K4 );
  for ( std::vector<Cell4>::const_iterator it = bdry.begin(), itE = bdry.end(); it != itE; ++it )
    {
      Cell4 cell = *it;
      Dimension d = K4.uDim( cell );
      CellMapConstIterator cmIt = complex4.findCell( d, cell );
      bdry_complex4.insertCell( d, cell, cmIt->second );
    }
  //bdry_complex4.insertCells( bdry.begin(), bdry.end() );
  trace.info() << "- [before collapse] K_bdry =" << bdry_complex4 << endl;
  functions::collapse( bdry_complex4, bdry.begin(), bdry.end(), priority, true, true, true );
  trace.info() << "- [after collapse]  K_bdry =" << bdry_complex4 << endl;
  for ( std::vector<Cell4>::const_iterator it = bdry.begin(), itE = bdry.end(); it != itE; ++it )
    {
      Cell4 cell  = *it;
      Dimension d = K4.uDim( cell );
      CellMapConstIterator cmIt = bdry_complex4.findCell( d, cell );
      if ( cmIt != bdry_complex4.end( d ) ) {
        CellMapIterator cmIt2 = complex4.findCell( d, cell );
        cmIt2->second = sure_data;
      }
    }
  trace.endBlock();

  //-------------- Collapse all -------------------------------------------
  trace.beginBlock( "Collapse all. " );
  std::copy( bdry.begin(), bdry.end(), std::back_inserter( inner ) );
  //typename CC4::DefaultCellMapIteratorPriority priority;
  functions::collapse( complex4, inner.begin(), inner.end(), priority, true, true, true );
  trace.info() << "- K=" << complex4 << endl;
  trace.endBlock();

  //-------------- Project complex onto surface --------------------------------
  trace.beginBlock( "Project complex onto surface. " );
  std::string project   = vm[ "project" ].as<std::string>();
  double epsilon        = vm[ "epsilon" ].as<double>();
  unsigned int max_iter = vm[ "max_iter" ].as<unsigned int>();
  std::vector<RealPoint3> points;
  if ( project == "Newton" )
    projectComplex( points, complex4, *shape, dshape, shape3, epsilon, max_iter );
  else
    doNotProjectComplex( points, complex4, *shape, dshape, shape3 );
  trace.endBlock();

  //-------------- Create Mesh -------------------------------------------
  trace.beginBlock( "Create Mesh. " );
  std::string view = vm[ "view" ].as<std::string>();
  bool highlight = ( view == "Singular" );
  Mesh<RealPoint3> mesh( true );
  std::map<Cell4,unsigned int> indices;
  int idx = 0;
  for ( CellMapConstIterator it = complex4.begin( 0 ), itE = complex4.end( 0 ); it != itE; ++it, ++idx )
    {
      Cell4 cell = it->first;
      indices[ cell ] = idx;
      mesh.addVertex( points[ idx ] );
    }
  for ( CellMapConstIterator it = complex4.begin( 2 ), itE = complex4.end( 2 ); it != itE; ++it, ++idx )
    {
      Cell4 cell = it->first;
      bool fixed = it->second.data & CC4::FIXED;
      Cells4 bdry = complex4.cellBoundary( cell, true );
      std::vector<unsigned int> face_idx;
      for ( Cells4::const_iterator itC = bdry.begin(), itCE = bdry.end(); itC != itCE; ++itC )
        {
          if ( complex4.dim( *itC ) == 0 )
            face_idx.push_back( indices[ *itC ] );
        }
      Color color = highlight
        ? ( fixed ? Color::White : Color(128,255,128)  )
        : Color::White;
      mesh.addQuadFace( face_idx[ 0 ], face_idx[ 1 ], face_idx[ 3 ], face_idx[ 2 ], color );
     }
  trace.endBlock();

  //-------------- View surface -------------------------------------------
  QApplication application(argc,argv);
  Point4 low4 = K4.lowerBound();
  Point4 up4  = K4.upperBound();
  KSpace3 K3;
  K3.init( Point3( low4[ 0 ], low4[ 1 ], low4[ 2 ] ),
           Point3( up4 [ 0 ], up4 [ 1 ], up4 [ 2 ] ), true );
  Viewer3D<Space3,KSpace3> viewer( K3 );
  viewer.setWindowTitle("Implicit surface viewer by 4d extension");
  viewer.show();
  viewer << mesh;
  viewer.setLineColor( highlight ? Color::Red : Color( 120, 120, 120 ) );
  // Drawing lines
  for ( CellMapConstIterator it = complex4.begin( 1 ), itE = complex4.end( 1 ); it != itE; ++it )
    {
      Cell4 cell = it->first;
      std::vector<Cell4> dummy;
      std::back_insert_iterator< std::vector<Cell4> > outIt1( dummy );
      complex4.directCoFaces( outIt1, cell );
      if ( ! dummy.empty() ) continue;
      Cells4 vertices = complex4.cellBoundary( cell );
      ASSERT( vertices.size() == 2 );
      RealPoint3 p1 = points[ indices[ vertices.front() ] ];
      RealPoint3 p2 = points[ indices[ vertices.back()  ] ];
      viewer.addLine( p1, p2, 0.05 );
    }
  viewer << Viewer3D<Space3,KSpace3>::updateDisplay;
  return application.exec();
}
///////////////////////////////////////////////////////////////////////////////
