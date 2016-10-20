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
#include "DGtal/images/CImage.h"
#include "DGtal/topology/CCellularGridSpaceND.h"
#include "DGtal/math/linalg/CLinearAlgebra.h"
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
      *
      * @tparam Calculus any discrete exterior calculus.
      * @tparam dim the dimension of the form.
      * @tparam duality either PRIMAL for a primal form or DUAL for a dual form.
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
      *
      * @tparam Calculus any discrete exterior calculus.
      * @tparam dim the dimension of the form.
      * @tparam duality either PRIMAL for a primal form or DUAL for a dual form.
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
      *
      * @tparam Calculus any discrete exterior calculus.
      * @tparam dim the dimension of the form.
      * @tparam duality either PRIMAL for a primal form or DUAL for a dual form.
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
      * @param[in,out] image any image of sufficient size.
      * @param pt a pixel coordinate (which is multiplied by \a pixel_size within).
      * @param val the value to write in \a pixel_size x \a pixel_size pixels.
      * @param pixel_size the chosen pixel_size (when 1, this is the normal setValue of an image).
      *
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Image>
      void writePixel( Image& image, typename Image::Point pt, typename Image::Value val,
                       int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      * @param[in,out] image any image of sufficient size.
      * @param pt a linel Khalimsky coordinates.
      * @param val the value to write in \a pixel_size x 1 pixels (if horizontal) or 1 x \a pixel_size pixels (if vertical).
      * @param pixel_size the chosen pixel_size (when 1, this is the normal setValue of an image).
      *
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Image>
      void writePrimalLinel( Image& image, typename Image::Point pt, typename Image::Value val,
                       int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      * @param[in,out] image any image of sufficient size.
      * @param pt a linel Khalimsky coordinates.
      * @param val the value to write in \a pixel_size x 1 pixels (if horizontal) or 1 x \a pixel_size pixels (if vertical).
      * @param pixel_size the chosen pixel_size (when 1, this is the normal setValue of an image).
      *
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Image>
      void writeDualLinel( Image& image, typename Image::Point pt, typename Image::Value val,
                           int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      
      /**
      * Displays the 2-form \a u in the given \a image. Scalar values
      * of \a u are first cut up and low according to \a cut_low and
      * \a cut_up, and then rescaled according to max and min
      * value. Then these values are transformed to image values with
      * the function \a functor. They are written in the image as
      * "pixels" of size \a pixel_size x \a pixel_size.
      *
      * @param calculus the discrete exterior calculus containing the 2-form \a u.
      * @param u any primal or dual 2-form defined in \a calculus.
      * @param[in,out] image the image where \a u is written.
      * @param functor the function transforming scalar values to image values.
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 2-form is mapped into \a image as \a pixel_size x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam AnyForm2 either a primal 2-form type or a dual 2-form type of the given Calculus.
      * @tparam Image any image type (see concepts::CImage).
      * @tparam Function any function type (double) -> typename Image::Value to convert form value to Image value.
      */
      template <typename Calculus, typename AnyForm2, typename Image, typename Function>
      void form2ToImage
      ( const Calculus& calculus, 
        const AnyForm2& u, 
        Image& image,
        const Function& functor,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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

      /**
      * Displays the primal or dual 1-form \a v in the given \a
      * image. Scalar values of \a v are first cut up and low
      * according to \a cut_low and \a cut_up, and then rescaled
      * according to max and min value. Then these values are
      * transformed to image values with the function \a functor. They
      * are written in the image as "lines" of size \a pixel_size x \a
      * 1 or \a 1 x \a pixel_size, depending on position and duality.
      *
      * @param calculus the discrete exterior calculus containing the 1-form \a v.
      * @param v any primal 1-form defined in \a calculus if \a primal is true, otherwise a dual 1-form.
      * @param primal tells if \a v is a primal 1-form (true), or a dual 1-form (false).
      * @param[in,out] image the image where \a v is written.
      * @param functor the function transforming scalar values to image values.
      * @param predicate the predicate telling for a value if it must be displayed (returns true in this case).
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 1-form is mapped into \a image as \a pixel_size x \a 1 pixels or \a 1 x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam Form1 either a primal 1-form if primal is true, or dual 1-form is primal is false.
      * @tparam Image any image type (see concepts::CImage).
      * @tparam Function any function type (double) -> typename Image::Value to convert form value to Image value.
      * @tparam Predicate any function type (double) -> bool to select 1-forms to display.
      */
      template <typename Calculus, typename Form1, typename Image, typename Function, typename Predicate>
      void form1ToImage
      ( const Calculus& calculus, 
        const Form1& v, bool primal,
        Image& image,
        const Function& functor,
        const Predicate& predicate,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      
      /**
      * Displays the dual 1-form \a v in the given \a image. Scalar
      * values of \a v are first cut up and low according to \a
      * cut_low and \a cut_up, and then rescaled according to max and
      * min value. Then these values are transformed to image values
      * with the function \a functor. They are written in the image as
      * "lines" of size \a pixel_size x \a 1 or \a 1 x \a pixel_size,
      * depending on position and duality.
      *
      * @param calculus the discrete exterior calculus containing the 1-form \a v.
      * @param v any dual 1-form defined in \a calculus.
      * @param[in,out] image the image where \a v is written.
      * @param functor the function transforming scalar values to image values.
      * @param predicate the predicate telling for a value if it must be displayed (returns true in this case).
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 1-form is mapped into \a image as \a pixel_size x \a 1 pixels or \a 1 x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam Image any image type (see concepts::CImage).
      * @tparam Function any function type (double) -> typename Image::Value to convert form value to Image value.
      * @tparam Predicate any function type (double) -> bool to select 1-forms to display.
      */
      template <typename Calculus, typename Image, typename Function, typename Predicate>
      void dualForm1ToImage
      ( const Calculus& calculus, 
        const typename Calculus::DualForm1& v, 
        Image& image,
        const Function& functor,
        const Predicate& predicate,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
        form1ToImage( calculus, v, false, image, functor, predicate,
                      cut_low, cut_up, pixel_size );
      }
      
      /**
      * Displays the primal 1-form \a v in the given \a image. Scalar
      * values of \a v are first cut up and low according to \a
      * cut_low and \a cut_up, and then rescaled according to max and
      * min value. Then these values are transformed to image values
      * with the function \a functor. They are written in the image as
      * "lines" of size \a pixel_size x \a 1 or \a 1 x \a pixel_size,
      * depending on position and duality.
      *
      * @param calculus the discrete exterior calculus containing the 1-form \a v.
      * @param v any primal 1-form defined in \a calculus.
      * @param[in,out] image the image where \a v is written.
      * @param functor the function transforming scalar values to image values.
      * @param predicate the predicate telling for a value if it must be displayed (returns true in this case).
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 1-form is mapped into \a image as \a pixel_size x \a 1 pixels or \a 1 x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam Image any image type (see concepts::CImage).
      * @tparam Function any function type (double) -> typename Image::Value to convert form value to Image value.
      * @tparam Predicate any function type (double) -> bool to select 1-forms to display.
      */
      template <typename Calculus, typename Image, typename Function, typename Predicate>
      void primalForm1ToImage
      ( const Calculus& calculus, 
        const typename Calculus::PrimalForm1& v, 
        Image& image,
        const Function& functor,
        const Predicate& predicate,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
        form1ToImage( calculus, v, true, image, functor, predicate,
                      cut_low, cut_up, pixel_size );
      }

      /**
      * Displays the three 2-forms \a u0, \a u1, \a u2 in the given \a image. Scalar values
      * of \a u0, \a u1, \a u2 are first cut up and low according to \a cut_low and
      * \a cut_up, and then rescaled according to max and min
      * value. Then these values are transformed to image values with
      * the function \a functor. They are written in the image as
      * "pixels" of size \a pixel_size x \a pixel_size.
      *
      * @param calculus the discrete exterior calculus containing the 2-forms \a u0, \a u1, \a u2.
      * @param u0 any primal or dual 2-form defined in \a calculus.
      * @param u1 any primal or dual 2-form defined in \a calculus.
      * @param u2 any primal or dual 2-form defined in \a calculus.
      * @param[in,out] image the image where \a u is written.
      * @param functor the function transforming three scalar values to image values.
      * @param cut_low every value of \a u0, \a u1, \a u2 below is set to \a cut_low.
      * @param cut_up  every value of \a u0, \a u1, \a u2 above is set to \a cut_up.
      * @param pixel_size every value of 2-forms is mapped into \a image as \a pixel_size x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam AnyForm2 either a primal 2-form type or a dual 2-form type of the given Calculus.
      * @tparam Image any image type (see concepts::CImage).
      * @tparam Function any function type (double,double,double) -> typename Image::Value to convert form value to Image value.
      */
      template <typename Calculus, typename AnyForm2, typename Image, typename Function>
      void threeForms2ToImage
      ( const Calculus& calculus, 
        const AnyForm2& u0, 
        const AnyForm2& u1, 
        const AnyForm2& u2, 
        Image& image,
        const Function& functor,
        double cut_low = 0.0, double cut_up = 1.0, int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      *
      * Displays the 2-form \a u in the given \a image. Scalar values
      * of \a u are first cut up and low according to \a cut_low and
      * \a cut_up, and then rescaled according to max and min
      * value. They are written in the image as
      * grey-level "pixels" of size \a pixel_size x \a pixel_size.
      *
      * @param calculus the discrete exterior calculus containing the 2-form \a u.
      * @param u any primal or dual 2-form defined in \a calculus.
      * @param[in,out] image the image where \a u is written.
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 2-form is mapped into \a image as \a pixel_size x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam AnyForm2 either a primal 2-form type or a dual 2-form type of the given Calculus.
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Calculus, typename AnyForm2, typename Image>
      void form2ToGreyLevelImage
      ( const Calculus& calculus, 
        const AnyForm2& u, 
        Image& image,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
        form2ToImage( calculus, u, image,
                      [] ( double x ) { return (unsigned char) ( round( x * 255.0 ) ); },
                      cut_low, cut_up, pixel_size );
      }

      /**
      * Standard method to output a primal 1-form into a grey-level image.
      *
      * Displays the primal 1-form \a v in the given \a image. Scalar
      * values of \a v are first cut up and low according to \a
      * cut_low and \a cut_up, and then rescaled according to max and
      * min value. If this value is belows 0.25, it is written in the image as
      * a "line" of size \a pixel_size x \a 1 or \a 1 x \a pixel_size,
      * depending on position.
      *
      * @param calculus the discrete exterior calculus containing the 1-form \a v.
      * @param v any primal 1-form defined in \a calculus.
      * @param[in,out] image the image where \a v is written.
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 1-form is mapped into \a image as \a pixel_size x \a 1 pixels or \a 1 x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Calculus, typename Image>
      void primalForm1ToGreyLevelImage
      ( const Calculus& calculus, 
        const typename Calculus::PrimalForm1& v, 
        Image& image,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      *
      * Displays the dual 1-form \a v in the given \a image. Scalar
      * values of \a v are first cut up and low according to \a
      * cut_low and \a cut_up, and then rescaled according to max and
      * min value. If this value is belows 0.25, it is written in the
      * image as a "line" of size \a pixel_size x \a 1 or \a 1 x \a
      * pixel_size, depending on position.
      *
      * @param calculus the discrete exterior calculus containing the 1-form \a v.
      * @param v any dual 1-form defined in \a calculus.
      * @param[in,out] image the image where \a v is written.
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 1-form is mapped into \a image as \a pixel_size x \a 1 pixels or \a 1 x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Calculus, typename Image>
      void dualForm1ToGreyLevelImage
      ( const Calculus& calculus, 
        const typename Calculus::DualForm1& v, 
        Image& image,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      *
      * Displays the primal 1-form \a v in the given \a image. Scalar
      * values of \a v are first cut up and low according to \a
      * cut_low and \a cut_up, and then rescaled according to max and
      * min value. If this value is belows 0.25, it is written in the image as
      * a "line" of size \a pixel_size x \a 1 or \a 1 x \a pixel_size,
      * depending on position, and of color \a color.
      *
      * @param calculus the discrete exterior calculus containing the 1-form \a v.
      * @param v any primal 1-form defined in \a calculus.
      * @param[in,out] image the image where \a v is written.
      * @param color the color for displaying 1-forms below 0.25.
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 1-form is mapped into \a image as \a pixel_size x \a 1 pixels or \a 1 x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Calculus, typename Image>
      void primalForm1ToRGBColorImage
      ( const Calculus& calculus, 
        const typename Calculus::PrimalForm1& v, 
        Image& image, Color color,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      *
      * Displays the dual 1-form \a v in the given \a image. Scalar
      * values of \a v are first cut up and low according to \a
      * cut_low and \a cut_up, and then rescaled according to max and
      * min value. If this value is belows 0.25, it is written in the image as
      * a "line" of size \a pixel_size x \a 1 or \a 1 x \a pixel_size,
      * depending on position, and of color \a color.
      *
      * @param calculus the discrete exterior calculus containing the 1-form \a v.
      * @param v any dual 1-form defined in \a calculus.
      * @param[in,out] image the image where \a v is written.
      * @param color the color for displaying 1-forms below 0.25.
      * @param cut_low every value of \a u below is set to \a cut_low.
      * @param cut_up  every value of \a u above is set to \a cut_up.
      * @param pixel_size every value of a 1-form is mapped into \a image as \a pixel_size x \a 1 pixels or \a 1 x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam Image any image type (see concepts::CImage).
      */
      template <typename Calculus, typename Image>
      void dualForm1ToRGBColorImage
      ( const Calculus& calculus, 
        const typename Calculus::DualForm1& v, 
        Image& image, Color color,
        double cut_low = 0.0, double cut_up = 1.0,
        int pixel_size = 1 )
      {
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
      *
      * Displays the three 2-forms \a u0, \a u1, \a u2 in the given \a
      * image as RGB colors. Scalar values of \a u0, \a u1, \a u2 are first cut up
      * and low according to \a cut_low and \a cut_up, and then
      * rescaled according to max and min value. Then these values are
      * transformed to RGB color image values (\a u0 defines the
      * intensity of the red channel, \a u1 the green channel, \a u2,
      * the blue channel). They are written in the image as "pixels"
      * of size \a pixel_size x \a pixel_size.
      *
      * @param calculus the discrete exterior calculus containing the 2-forms \a u0, \a u1, \a u2.
      * @param u0 any primal or dual 2-form defined in \a calculus.
      * @param u1 any primal or dual 2-form defined in \a calculus.
      * @param u2 any primal or dual 2-form defined in \a calculus.
      * @param[in,out] image the image where \a u is written.
      * @param cut_low every value of \a u0, \a u1, \a u2 below is set to \a cut_low.
      * @param cut_up  every value of \a u0, \a u1, \a u2 above is set to \a cut_up.
      * @param pixel_size every value of 2-forms is mapped into \a image as \a pixel_size x \a pixel_size pixels.
      *
      * @tparam Calculus any discrete exterior calculus type.
      * @tparam AnyForm2 either a primal 2-form type or a dual 2-form type of the given Calculus.
      * @tparam Image any image type (see concepts::CImage).
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
        BOOST_CONCEPT_ASSERT(( concepts::CImage<Image> ));
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
  * Description of template class 'DECImage2D' <p> \brief Aim: This
  * class simplifies the development of 2D image processing tools
  * using discrete exterior calculus. Most notably it take care of
  * initializing correctly a discrete exterior calculus in some 2D
  * domain, and precomputes derivative and Hodge star operators.  You
  * may have a look at module \ref moduleAT for such image processing
  * tools and more explanation on discrete calculus.
  *
  * @tparam TKSpace any model of CCellularGridSpaceND, e.g KhalimskySpaceND
  * @tparam TLinearAlgebra any back-end for performing linear algebra, default is EigenLinearAlgebraBackend.
  *
  * @see ATu0v1
  * @see ATu2v0
  *
  */
  template < typename TKSpace,
             typename TLinearAlgebra = EigenLinearAlgebraBackend >
  struct DECImage2D {
    typedef TKSpace                                        KSpace;
    typedef TLinearAlgebra                                 LinearAlgebra;
    BOOST_CONCEPT_ASSERT(( concepts::CCellularGridSpaceND< KSpace > )); 
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
    * Constructor from cellular grid space, which specifies the domain of calculus.
    *
    * @param aKSpace the cellular grid space specifies the domain of
    * calculus (i.e. all the cells and incidence), which is cloned
    * inside the class.
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
