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

#pragma once

/**
 * @file
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Marion Foare (\c marion.foare@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2016/10/12
 *
 * Header file for module DECImageHelpers.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DECImageHelpers_RECURSES)
#error Recursive header files inclusion detected in DECImageHelpers.h
#else // defined(DECImageHelpers_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DECImageHelpers_RECURSES

#if !defined DECImageHelpers_h
/** Prevents repeated inclusion of headers. */
#define DECImageHelpers_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/dec/DiscreteExteriorCalculus.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal {
  namespace functions {
    namespace dec {
      
      /**
      * Builds a diagonal linear operator from a k-form. These
      * operators arise naturally when differentiating with respect to
      * another variable (e.g. d/dx (vx)^t (vx) = diag(v^2) x).
      *
      * @param[in] kform any kform w.
      * @return the corresponding linear operator diag(w)
      */
      template <typename Calculus, DGtal::Dimension dim, DGtal::Duality duality> 
      DGtal::LinearOperator<Calculus, dim, duality, dim, duality> 
      diagonal(const DGtal::KForm<Calculus, dim, duality>& kform) 
      { 
        typedef DGtal::LinearOperator<Calculus,dim, duality, dim, duality> Operator; 
        typedef typename Calculus::LinearAlgebraBackend::Triplet           Triplet; 
        typedef typename Calculus::Index                                   Index; 
        typedef std::vector<Triplet>                                       Triplets;

        Triplets triplets;
        for (Index index=0; index<kform.length(); index++) 
          triplets.push_back(Triplet(index, index, kform.myContainer(index)));
        
        Operator op( kform.myCalculus );
        op.myContainer.setFromTriplets( triplets.begin(), triplets.end() );
        return op;
      }

      /**
      * Squares the given k-form.
      *
      * @param[in,out] kform any kform.
      */
      template <typename Calculus, DGtal::Dimension dim, DGtal::Duality duality> 
      void
      squares(DGtal::KForm<Calculus, dim, duality>& kform) 
      {
        kform.myContainer.array() = kform.myContainer.array().square();
      }
      
      /**
      * Builds a diagonal linear operator from a k-form and squares it. These
      * operators arise naturally when differentiating with respect to
      * another variable (e.g. d/dx (vx)^t (vx) = diag(v^2) x).
      *
      * @param[in] kform any kform v
      * @return the corresponding linear operator diag(v^2)
      */
      template <typename Calculus, DGtal::Dimension dim, DGtal::Duality duality> 
      DGtal::LinearOperator<Calculus, dim, duality, dim, duality> 
      squaredDiagonal(const DGtal::KForm<Calculus, dim, duality>& kform) 
      {
        auto v2 = kform;
        squares( v2 );
        return diagonal( v2 );
      }

      /**
      * Considers an image \a image to have pixels of size \a
      * pixel_size x \a pixel_size, and writes the value \a val at the
      * specified pixel position \a pt.
      *
      * @param[in,out] any image of sufficient size.
      * @param pt a pixel coordinate (which is multiplied by \a pixel_size within).
      * @param val the value to write in \a pixel_size x \a pixel_size pixels.
      * @param pixel_size the chosen pixel_size (when 1, this is the normal setValue of an image).
      */
      template <typename Image>
      void writePixel( Image& image, typename Image::Point pt, typename Image::Value val,
                       int pixel_size = 1 )
      {
        typedef typename Image::Point      Point;
        typedef typename Point::Coordinate Coordinate;
        pt *= pixel_size;
        for ( int y = 0; y < pixel_size; y++ )
          for ( int x = 0; x < pixel_size; x++ )
            {
              Point q( (Coordinate) x, (Coordinate) y );
              image.setValue( pt + q, val );
            }
      }

      /**
      * Considers an image \a image to have pixels of size \a
      * pixel_size x \a pixel_size, and writes the value \a val at the
      * specified linel position \a pt.
      *
      * @param[in,out] any image of sufficient size.
      * @param pt a linel Khalimsky coordinates.
      * @param val the value to write in \a pixel_size x 1 pixels (if horizontal) or 1 x \a pixel_size pixels (if vertical).
      * @param pixel_size the chosen pixel_size (when 1, this is the normal setValue of an image).
      */
      template <typename Image>
      void writePrimalLinel( Image& image, typename Image::Point pt, typename Image::Value val,
                       int pixel_size = 1 )
      {
        typedef typename Image::Point      Point;
        typedef typename Point::Coordinate Coordinate;
        int pixel_size_x = NumberTraits<Coordinate>::even( pt[ 0 ] ) ? 1 : pixel_size;
        int pixel_size_y = NumberTraits<Coordinate>::even( pt[ 1 ] ) ? 1 : pixel_size;
        pt /= 2;
        pt *= pixel_size;
        for ( int y = 0; y < pixel_size_y; y++ )
          for ( int x = 0; x < pixel_size_x; x++ )
            {
              Point q( (Coordinate) x, (Coordinate) y );
              image.setValue( pt + q, val );
            }
      }
      
      /**
      * Considers an image \a image to have pixels of size \a
      * pixel_size x \a pixel_size, and writes the value \a val at the
      * specified linel position \a pt.
      *
      * @param[in,out] any image of sufficient size.
      * @param pt a linel Khalimsky coordinates.
      * @param val the value to write in \a pixel_size x 1 pixels (if horizontal) or 1 x \a pixel_size pixels (if vertical).
      * @param pixel_size the chosen pixel_size (when 1, this is the normal setValue of an image).
      */
      template <typename Image>
      void writeDualLinel( Image& image, typename Image::Point pt, typename Image::Value val,
                           int pixel_size = 1 )
      {
        typedef typename Image::Point      Point;
        typedef typename Point::Coordinate Coordinate;
        int pixel_size_x = NumberTraits<Coordinate>::even( pt[ 0 ] ) ? 0 : pixel_size-1;
        int pixel_size_y = NumberTraits<Coordinate>::even( pt[ 1 ] ) ? 0 : pixel_size-1;
        pt /= 2;
        pt *= pixel_size;
        for ( int y = pixel_size_y; y < pixel_size; y++ )
          for ( int x = pixel_size_x; x < pixel_size; x++ )
            {
              Point q( (Coordinate) x, (Coordinate) y );
              image.setValue( pt + q, val );
            }
      }
      
      template <typename Calculus, typename AnyForm2, typename Image>
      void form2ToImage
      ( const Calculus& calculus, 
        const AnyForm2& u, 
        Image& image,
        std::function< typename Image::Value( double ) > functor,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        typedef typename Calculus::Index  Index;
        typedef typename Calculus::SCell  SCell;
        typedef typename Calculus::Scalar Scalar;
        typedef typename Calculus::KSpace KSpace;
        typedef typename KSpace::Point    Point;
        typedef typename KSpace::Integer  Integer;
        double min_u = NumberTraits<Scalar>::castToDouble( u.myContainer[ 0 ] );
        double max_u = min_u;
        for ( Index index = 0; index < u.myContainer.rows(); index++)
          {
            double v = NumberTraits<Scalar>::castToDouble( u.myContainer[ index ] );
            min_u = std::min( min_u, v );
            max_u = std::max( max_u, v );
          }
        if ( min_u < cut_low ) min_u = cut_low;
        if ( max_u > cut_up  ) max_u = cut_up;
        for ( Index index = 0; index < u.myContainer.rows(); index++)
          {
            SCell cell = u.getSCell( index );
            double v = NumberTraits<Scalar>::castToDouble( u.myContainer[ index ] );
            double w = std::min( cut_up, std::max( cut_low, v ) );
            if ( min_u != max_u ) w = ( w - min_u ) / ( max_u - min_u );
            writePixel( image, calculus.myKSpace.sCoords( cell ), functor( w ), pixel_size );
          }
      }

      // template <typename Calculus, typename Image>
      // void dualForm2ToImage
      // ( const Calculus& calculus, 
      //   const typename Calculus::DualForm2& u, 
      //   Image& image,
      //   std::function< typename Image::Value( double ) > functor,
      //   double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      // {
      //   typedef typename Calculus::Index  Index;
      //   typedef typename Calculus::SCell  SCell;
      //   typedef typename Calculus::Scalar Scalar;
      //   typedef typename Calculus::KSpace KSpace;
      //   typedef typename KSpace::Point    Point;
      //   typedef typename KSpace::Integer  Integer;
      //   double min_u = NumberTraits<Scalar>::castToDouble( u.myContainer[ 0 ] );
      //   double max_u = min_u;
      //   for ( Index index = 0; index < u.myContainer.rows(); index++)
      //     {
      //       double v = NumberTraits<Scalar>::castToDouble( u.myContainer[ index ] );
      //       min_u = std::min( min_u, v );
      //       max_u = std::max( max_u, v );
      //     }
      //   if ( min_u < cut_low ) min_u = cut_low;
      //   if ( max_u > cut_up  ) max_u = cut_up;
      //   for ( Index index = 0; index < u.myContainer.rows(); index++)
      //     {
      //       SCell cell = u.getSCell( index );
      //       double v = NumberTraits<Scalar>::castToDouble( u.myContainer[ index ] );
      //       double w = std::min( cut_up, std::max( cut_low, v ) );
      //       if ( min_u != max_u ) w = ( w - min_u ) / ( max_u - min_u );
      //       writePixel( image, calculus.myKSpace.sCoords( cell ), functor( w ), pixel_size );
      //     }
      // }

      template <typename Calculus, typename Form1, typename Image>
      void form1ToImage
      ( const Calculus& calculus, 
        const Form1& v, bool primal,
        Image& image,
        std::function< typename Image::Value( double ) > functor,
        std::function< bool ( double ) > predicate,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        typedef typename Calculus::Index  Index;
        typedef typename Calculus::SCell  SCell;
        typedef typename Calculus::Scalar Scalar;
        typedef typename Calculus::KSpace KSpace;
        typedef typename KSpace::Point    Point;
        typedef typename KSpace::Integer  Integer;
        double min_v = NumberTraits<Scalar>::castToDouble( v.myContainer[ 0 ] );
        double max_v = min_v;
        for ( Index index = 0; index < v.myContainer.rows(); index++)
          {
            double w = NumberTraits<Scalar>::castToDouble( v.myContainer[ index ] );
            min_v = std::min( min_v, w );
            max_v = std::max( max_v, w );
          }
        if ( min_v < cut_low ) min_v = cut_low;
        if ( max_v > cut_up  ) max_v = cut_up;
        for ( Index index = 0; index < v.myContainer.rows(); index++)
          {
            SCell cell = v.getSCell( index );
            double u = NumberTraits<Scalar>::castToDouble( v.myContainer[ index ] );
            if ( ! predicate( u ) ) continue; 
            double w = std::min( cut_up, std::max( cut_low, u ) );
            if ( min_v != max_v ) w = ( w - min_v ) / ( max_v - min_v );
            Point kpt = calculus.myKSpace.sKCoords( cell );
            if ( primal ) writePrimalLinel( image, kpt, functor( w ), pixel_size );
            else          writeDualLinel  ( image, kpt, functor( w ), pixel_size );
          }
      }
      
      template <typename Calculus, typename Image>
      void dualForm1ToImage
      ( const Calculus& calculus, 
        const typename Calculus::DualForm1& v, 
        Image& image,
        std::function< typename Image::Value( double ) > functor,
        std::function< bool ( double ) > predicate,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        form1ToImage( calculus, v, false, image, functor, predicate,
                      cut_low, cut_up, pixel_size );
      }
      
      template <typename Calculus, typename Image>
      void primalForm1ToImage
      ( const Calculus& calculus, 
        const typename Calculus::PrimalForm1& v, 
        Image& image,
        std::function< typename Image::Value( double ) > functor,
        std::function< bool ( double ) > predicate,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        form1ToImage( calculus, v, true, image, functor, predicate,
                      cut_low, cut_up, pixel_size );
      }

      template <typename Calculus, typename AnyForm2, typename Image>
      void threeForms2ToImage
      ( const Calculus& calculus, 
        const AnyForm2& u0, 
        const AnyForm2& u1, 
        const AnyForm2& u2, 
        Image& image,
        std::function< typename Image::Value( double, double, double ) > functor,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        typedef typename Calculus::Index  Index;
        typedef typename Calculus::SCell  SCell;
        typedef typename Calculus::Scalar Scalar;
        double min_u = NumberTraits<Scalar>::castToDouble( u0.myContainer[ 0 ] );
        double max_u = min_u;
        for ( Index index = 0; index < u0.myContainer.rows(); index++)
          {
            double v = NumberTraits<Scalar>::castToDouble( u0.myContainer[ index ] );
            min_u = std::min( min_u, v );
            max_u = std::max( max_u, v );
          }
        for ( Index index = 0; index < u1.myContainer.rows(); index++)
          {
            double v = NumberTraits<Scalar>::castToDouble( u1.myContainer[ index ] );
            min_u = std::min( min_u, v );
            max_u = std::max( max_u, v );
          }
        for ( Index index = 0; index < u2.myContainer.rows(); index++)
          {
            double v = NumberTraits<Scalar>::castToDouble( u2.myContainer[ index ] );
            min_u = std::min( min_u, v );
            max_u = std::max( max_u, v );
          }
        if ( min_u < cut_low ) min_u = cut_low;
        if ( max_u > cut_up  ) max_u = cut_up;
        for ( Index index = 0; index < u0.myContainer.rows(); index++)
          {
            SCell cell = u0.getSCell( index );
            double v0 = NumberTraits<Scalar>::castToDouble( u0.myContainer[ index ] );
            double w0 = std::min( cut_up, std::max( cut_low, v0 ) );
            if ( min_u != max_u ) w0 = ( w0 - min_u ) / ( max_u - min_u );
            double v1 = NumberTraits<Scalar>::castToDouble( u1.myContainer[ index ] );
            double w1 = std::min( cut_up, std::max( cut_low, v1 ) );
            if ( min_u != max_u ) w1 = ( w1 - min_u ) / ( max_u - min_u );
            double v2 = NumberTraits<Scalar>::castToDouble( u2.myContainer[ index ] );
            double w2 = std::min( cut_up, std::max( cut_low, v2 ) );
            if ( min_u != max_u ) w2 = ( w2 - min_u ) / ( max_u - min_u );
            writePixel( image, calculus.myKSpace.sCoords( cell ), functor( w0, w1, w2 ), pixel_size );
          }
      }
      
      /**
      * Standard method to output a 2-form into a grey-level image.
      */
      template <typename Calculus, typename AnyForm2, typename Image>
      void form2ToGreyLevelImage
      ( const Calculus& calculus, 
        const AnyForm2& u, 
        Image& image,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        form2ToImage( calculus, u, image,
                      [] ( double x ) { return (unsigned char) ( round( x * 255.0 ) ); },
                      cut_low, cut_up, pixel_size );
      }

      /**
      * Standard method to output a primal 1-form into a grey-level image.
      */
      template <typename Calculus, typename Image>
      void primalForm1ToGreyLevelImage
      ( const Calculus& calculus, 
        const typename Calculus::PrimalForm1& v, 
        Image& image,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        // Threshold is 0.25 instead of 0.5 because an edge connecting
        // two vertices with v=0 and v=1 should not belong to the
        // discontinuity set.
        primalForm1ToImage( calculus, v, image,
                            [] ( double x ) { return (unsigned char) ( round( x * 255.0 ) ); },
                            [] ( double x ) { return x < 0.25; },
                            cut_low, cut_up, pixel_size );
      }
      
      /**
      * Standard method to output a dual 1-form into a grey-level image.
      */
      template <typename Calculus, typename Image>
      void dualForm1ToGreyLevelImage
      ( const Calculus& calculus, 
        const typename Calculus::DualForm1& v, 
        Image& image,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        // Threshold is 0.25 instead of 0.5 because an edge connecting
        // two vertices with v=0 and v=1 should not belong to the
        // discontinuity set.
        dualForm1ToImage( calculus, v, image,
                          [] ( double x ) { return (unsigned char) ( round( x * 255.0 ) ); },
                          [] ( double x ) { return x < 0.25; },
                          cut_low, cut_up, pixel_size );
      }

      /**
      * Standard method to output a primal 1-form into a color image.
      */
      template <typename Calculus, typename Image>
      void primalForm1ToRGBColorImage
      ( const Calculus& calculus, 
        const typename Calculus::PrimalForm1& v, 
        Image& image, Color color,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        // Threshold is 0.25 instead of 0.5 because an edge connecting
        // two vertices with v=0 and v=1 should not belong to the
        // discontinuity set.
        primalForm1ToImage( calculus, v, image,
                            [color] ( double x ) { return color; },
                            [] ( double x ) { return x < 0.25; },
                            cut_low, cut_up, pixel_size );
      }

      /**
      * Standard method to output a dual 1-form into a color image.
      */
      template <typename Calculus, typename Image>
      void dualForm1ToRGBColorImage
      ( const Calculus& calculus, 
        const typename Calculus::DualForm1& v, 
        Image& image, Color color,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        // Threshold is 0.25 instead of 0.5 because an edge connecting
        // two vertices with v=0 and v=1 should not belong to the
        // discontinuity set.
        dualForm1ToImage( calculus, v, image,
                          [color] ( double x ) { return color; },
                          [] ( double x ) { return x < 0.25; },
                          cut_low, cut_up, pixel_size );
      }

      /**
      * Standard method to output three 2-forms into a RGB Color image.
      */
      template <typename Calculus, typename AnyForm2, typename Image>
      void threeForms2ToRGBColorImage
      ( const Calculus& calculus, 
        const AnyForm2& u0, 
        const AnyForm2& u1, 
        const AnyForm2& u2, 
        Image& image,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        threeForms2ToImage
          ( calculus, u0, u1, u2, image,
            [] ( double r, double g, double b )
            { return Color( (unsigned char) ( round( r * 255.0 ) ),
                            (unsigned char) ( round( g * 255.0 ) ),
                            (unsigned char) ( round( b * 255.0 ) ) ); },
            cut_low, cut_up, pixel_size );
      }

    } // namespace dec
  } // namespace functions

  /////////////////////////////////////////////////////////////////////////////
  // template class DECImage2D
  /**
  * Description of template class 'DECImage2D' <p> \brief Aim: This class
  * simplifies the development of 2D image processing tools using discrete exterior calculus.
  *
  */
  template < typename TKSpace,
             typename TLinearAlgebra = EigenLinearAlgebraBackend >
  struct DECImage2D {
    typedef TKSpace                                        KSpace;
    typedef TLinearAlgebra                                 LinearAlgebra;
    typedef typename KSpace::Space                         Space;
    typedef typename Space::Point                          Point;
    typedef typename Space::RealVector                     RealVector;
    typedef typename RealVector::Component                 Scalar;
    typedef typename KSpace::SCell                         SCell;
    typedef typename KSpace::Cell                          Cell;
    typedef typename KSpace::Surfel                        Surfel;
    typedef HyperRectDomain<Space>                         Domain;
    typedef DiscreteExteriorCalculus<2,2, LinearAlgebra>   Calculus;
    typedef DiscreteExteriorCalculusFactory<LinearAlgebra> CalculusFactory;
    typedef typename Calculus::Index                       Index;
    typedef typename Calculus::PrimalForm0                 PrimalForm0;
    typedef typename Calculus::PrimalForm1                 PrimalForm1;
    typedef typename Calculus::PrimalForm2                 PrimalForm2;
    typedef typename Calculus::PrimalIdentity0             PrimalIdentity0;
    typedef typename Calculus::PrimalIdentity1             PrimalIdentity1;
    typedef typename Calculus::PrimalIdentity2             PrimalIdentity2;
    typedef typename Calculus::PrimalDerivative0           PrimalDerivative0;
    typedef typename Calculus::PrimalDerivative1           PrimalDerivative1;
    typedef typename Calculus::DualDerivative0             DualDerivative0;
    typedef typename Calculus::DualDerivative1             DualDerivative1;
    typedef typename Calculus::PrimalAntiderivative1       PrimalAntiderivative1;
    typedef typename Calculus::PrimalAntiderivative2       PrimalAntiderivative2;
    typedef typename Calculus::DualAntiderivative1         DualAntiderivative1;
    typedef typename Calculus::DualAntiderivative2         DualAntiderivative2;
    typedef typename Calculus::PrimalHodge0                PrimalHodge0;
    typedef typename Calculus::PrimalHodge1                PrimalHodge1;
    typedef typename Calculus::PrimalHodge2                PrimalHodge2;
    typedef typename Calculus::DualHodge0                  DualHodge0;
    typedef typename Calculus::DualHodge1                  DualHodge1;
    typedef typename Calculus::DualHodge2                  DualHodge2;
    typedef typename LinearAlgebra::SolverSimplicialLLT    LinearAlgebraSolver;
    typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 0, PRIMAL, 0, PRIMAL> 
                                                           SolverU;
    typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 1, PRIMAL, 1, PRIMAL> 
                                                           SolverV;

    BOOST_STATIC_ASSERT(( KSpace::dimension == 2 ));

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~DECImage2D() = default;

    /**
    * Default constructor. The object needs to be initialized with \ref init.
    * @param _verbose specifies the verbose level (0: silent, 1: more info ... ). 
    */
    DECImage2D( int _verbose = 1 )
      : verbose( _verbose ), 
        calculus(), 
        D0( calculus ), D1( calculus ), 
        dual_D0( calculus ), dual_D1( calculus ),
        primal_h0( calculus ), primal_h1( calculus ), primal_h2( calculus ),
        dual_h0( calculus ), dual_h1( calculus ), dual_h2( calculus )
    {}
    
    /**
    * Constructor from Khalimsky space, which specifies the domain of calculus.
    */
    void init( Clone<KSpace> aKSpace )
    {
      calculus.myKSpace = aKSpace;
      const KSpace & K  = calculus.myKSpace;
      domain            = Domain( K.lowerBound(), K.upperBound() );
      Point  p0         = K.uKCoords( K.lowerCell() );
      Point  p1         = K.uKCoords( K.upperCell() );
      cell_domain       = Domain( p0, p1 );

      if ( verbose > 0 ) trace.beginBlock("building AT functionnals");
      // Adds all the cell
      for ( typename Domain::ConstIterator it = cell_domain.begin(), itE = cell_domain.end(); 
            it != itE; ++it )
        calculus.insertSCell( K.sCell( *it ) ); // ajoute toutes les cellules de Khalimsky.
      calculus.updateIndexes();
      if ( verbose > 1 ) trace.info() << calculus << std::endl;
      // Precomputes operators.
      if ( verbose > 1 ) trace.info() << "primal_D0" << std::endl;
      D0 = calculus.template derivative<0,PRIMAL>();
      if ( verbose > 1 ) trace.info() << "primal_D1" << std::endl;
      D1 = calculus.template derivative<1,PRIMAL>();
      if ( verbose > 1 ) trace.info() << "dual_D0" << std::endl;
      dual_D0   = calculus.template derivative<0,DUAL>();
      if ( verbose > 1 ) trace.info() << "dual_D1" << std::endl;
      dual_D1   = calculus.template derivative<1,DUAL>();
      if ( verbose > 1 ) trace.info() << "primal_h1" << std::endl;
      primal_h0 = calculus.template hodge<0,PRIMAL>();
      if ( verbose > 1 ) trace.info() << "primal_h1" << std::endl;
      primal_h1 = calculus.template hodge<1,PRIMAL>();
      if ( verbose > 1 ) trace.info() << "primal_h2" << std::endl;
      primal_h2 = calculus.template hodge<2,PRIMAL>();
      if ( verbose > 1 ) trace.info() << "dual_h1" << std::endl;
      dual_h0   = calculus.template hodge<0,DUAL>();
      if ( verbose > 1 ) trace.info() << "dual_h1" << std::endl;
      dual_h1   = calculus.template hodge<1,DUAL>();
      if ( verbose > 1 ) trace.info() << "dual_h2" << std::endl;
      dual_h2   = calculus.template hodge<2,DUAL>();
      if ( verbose > 0 ) trace.endBlock();
    }

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    DECImage2D ( const DECImage2D & other ) = delete;

    /**
     * Move constructor.
     * @param other the object to move.
     */
    DECImage2D ( DECImage2D && other ) = delete;

    /**
     * Copy assignment operator.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    DECImage2D & operator= ( const DECImage2D & other ) = delete;

    /**
     * Move assignment operator.
     * @param other the object to move.
     * @return a reference on 'this'.
     */
    DECImage2D & operator= ( DECImage2D && other ) = delete;

    // ------------------------- Public Datas ------------------------------
  public:

    /// The verbose level (0: silent).
    int verbose;
    /// The discrete exterior calculus instance.
    Calculus calculus;
    /// The image domain (i.e. all the pixels)
    Domain   domain;
    /// The cell domain (i.e. all the cells)
    Domain   cell_domain;
    /// primal derivative: 0-form -> 1-form
    PrimalDerivative0 D0;
    /// primal derivative: 1-form -> 2-form
    PrimalDerivative1 D1;
    /// dual derivative dual 0-form -> dual 1-form
    DualDerivative0   dual_D0;
    /// dual derivative dual 1-form -> dual 2-form
    DualDerivative1   dual_D1;
    /// hodge star: 0-form -> dual 0-form
    PrimalHodge0      primal_h0;
    /// hodge star: 1-form -> dual 1-form
    PrimalHodge1      primal_h1;
    /// hodge star: 2-form -> dual 2-form
    PrimalHodge2      primal_h2;
    /// hodge star: dual 0-form -> 0-form
    DualHodge0        dual_h0;
    /// hodge star: dual 1-form -> 1-form
    DualHodge1        dual_h1;
    /// hodge star: dual 2-form -> 2-form
    DualHodge2        dual_h2;

  };

} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DECImageHelpers_h

#undef DECImageHelpers_RECURSES
#endif // else defined(DECImageHelpers_RECURSES)
