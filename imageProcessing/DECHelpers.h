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
 *
 * @date 2016/10/12
 *
 * Header file for module DECHelpers.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DECHelpers_RECURSES)
#error Recursive header files inclusion detected in DECHelpers.h
#else // defined(DECHelpers_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DECHelpers_RECURSES

#if !defined DECHelpers_h
/** Prevents repeated inclusion of headers. */
#define DECHelpers_h

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


    } // namespace dec
  } // namespace functions
} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DECHelpers_h

#undef DECHelpers_RECURSES
#endif // else defined(DECHelpers_RECURSES)
