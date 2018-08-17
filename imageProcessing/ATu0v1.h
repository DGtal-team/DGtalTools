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
 * @file ATu0v1.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Marion Foare (\c marion.foare@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2016/10/12
 *
 * Header file for module ATu0v1.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ATu0v1_RECURSES)
#error Recursive header files inclusion detected in ATu0v1.h
#else // defined(ATu0v1_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ATu0v1_RECURSES

#if !defined ATu0v1_h
/** Prevents repeated inclusion of headers. */
#define ATu0v1_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
// always include EigenSupport.h before any other Eigen headers
#include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/base/Common.h"
#include "DGtal/dec/DiscreteExteriorCalculus.h"
#include "DGtal/dec/DiscreteExteriorCalculusSolver.h"
#include "DGtal/dec/DiscreteExteriorCalculusFactory.h"
#include "DECImageHelpers.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class ATu0v1
  /**
  * Description of template class 'ATu0v1' <p> \brief Aim: This class
  * solves Ambrosio-Tortorelli functional in a plane for \a u a
  * (vector of) 0-form(s) and \a v a 1-form. \a u is a regularized
  * approximation of an input data \a g, while \a v represents the
  * set of discontinuities of \a u.
  *
  * @tparam TKSpace any model of CCellularGridSpaceND, e.g KhalimskySpaceND
  * @tparam TLinearAlgebra any back-end for performing linear algebra, default is EigenLinearAlgebraBackend.
  *
  */
  template < typename TKSpace,
             typename TLinearAlgebra = EigenLinearAlgebraBackend >
  struct ATu0v1 : public DECImage2D<TKSpace, TLinearAlgebra>
  {
    typedef TKSpace                             KSpace;
    typedef TLinearAlgebra                      LinearAlgebra;
    typedef DECImage2D<TKSpace, TLinearAlgebra> Base;
    using typename Base::Space;
    using typename Base::Point;
    using typename Base::RealVector;
    using typename Base::Scalar;
    using typename Base::SCell;
    using typename Base::Domain;
    using typename Base::Calculus;
    using typename Base::Index;
    using typename Base::PrimalForm0;
    using typename Base::PrimalForm1;
    using typename Base::PrimalForm2;
    using typename Base::PrimalIdentity0;
    using typename Base::PrimalIdentity1;
    using typename Base::PrimalIdentity2;
    using typename Base::PrimalDerivative0;
    using typename Base::PrimalDerivative1;
    using typename Base::DualDerivative0;
    using typename Base::DualDerivative1;
    using typename Base::PrimalHodge0;
    using typename Base::PrimalHodge1;
    using typename Base::PrimalHodge2;
    using typename Base::DualHodge0;
    using typename Base::DualHodge1;
    using typename Base::DualHodge2;
    using Base::calculus;
    using Base::verbose;
    using Base::D0;
    using Base::D1;
    using Base::dual_D0;
    using Base::dual_D1;
    using Base::primal_h0;
    using Base::primal_h1;
    using Base::primal_h2;
    using Base::dual_h0;
    using Base::dual_h1;
    using Base::dual_h2;

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
    ~ATu0v1() = default;

    /**
    * Default constructor. The object needs to be initialized with \ref init.
    * @param _verbose specifies the verbose level (0: silent, 1: more info ... ). 
    */
    ATu0v1( int _verbose = 1 );
    
    /**
    * Constructor from Khalimsky space, which specifies the domain of calculus.
    */
    void init( Clone<KSpace> K );

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    ATu0v1 ( const ATu0v1 & other ) = delete;

    /**
     * Move constructor.
     * @param other the object to move.
     */
    ATu0v1 ( ATu0v1 && other ) = delete;

    /**
     * Copy assignment operator.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    ATu0v1 & operator= ( const ATu0v1 & other ) = delete;

    /**
     * Move assignment operator.
     * @param other the object to move.
     * @return a reference on 'this'.
     */
    ATu0v1 & operator= ( ATu0v1 && other ) = delete;

    /**
    * Adds an input 0-form by filtering an \a image values.
    *
    * @param image any image such that the domain of this space is
    * included in the domain of the image.
    *
    * @param f any functor associating a scalar to an image value.
    *
    * @param perfect_data when 'false', this is normal input data,
    * otherwise this is perfect data only used for SNR computation.
    *
    * @note For a grey-level image stored with values `unsigned char`,
    * should be called as
    *
    * @code
    * AT.addInput( image, [] (unsigned char c ) { return (double) c / 255.0; } );
    * @endcode
    *
    * @note For a color image stored with values `Color`, should be called as
    * @code
    * AT.addInput( image, [] ( Color c ) { return (double) c.red()   / 255.0; } );
    * AT.addInput( image, [] ( Color c ) { return (double) c.green() / 255.0; } );
    * AT.addInput( image, [] ( Color c ) { return (double) c.blue()  / 255.0; } );
    * @endcode
    *
    * @tparam Image any Image type.
    * @tparam Function any function type ( typename Image::Value ) -> Scalar.
    */
    template <typename Image, typename Function >
    void addInput( const Image& image,
                   const Function& f,
                   bool perfect_data = false );

    /// Sets approximation \a u to be equal to the input. Used for
    /// initializating \a u. Should be called once all \ref addInput
    /// have been called.
    void setUFromInput();

    /** 
    * Sets the parameter \f$ alpha \f$ as global to the image. Should be
    * set \b before \ref setLambda and \ref setEpsilon.
    *
    * @param _alpha the \f$ \alpha \f$ parameter in AT functional ( in term \f$
    * \int \alpha | u - g |^2 \f$ ). Dimension theory tells that it is in
    * 1/area unit, the lower the smoother will be the output.
    */
    void setAlpha( Scalar _alpha );

    /** 
    * Sets the parameter \f$ alpha \f$ of the image, as well as a
    * weight for each input data as the primal 0-form m. Should be set
    * \b before \ref setLambda and \ref setEpsilon.
    *
    * @note Useful for inpainting applications where you indicate with
    * m=0 that the specified pixel data is useless.
    *
    * @param _alpha the \f$ \alpha \f$ parameter in AT functional ( in term \f$
    * \int \alpha | u - g |^2 \f$ ). Dimension theory tells that it is in
    * 1/area unit, the lower the smoother will be the output.
    *
    * @param m a 0-form that specifies which input data is significant
    * (1) or not be used (0).
    */
    void setAlpha( Scalar _alpha, const PrimalForm0& m );

    /**
    * Sets the parameter \f$ \lambda \f$ of AT functional. Should be
    * set \b after \ref setAlpha and \b before \ref setEpsilon.
    *
    * @param _lambda the \f$ \lambda \f$ parameter in AT functional (
    * in terms \f$ \int \lambda \epsilon | v \grad u |^2 + \int
    * \frac{\lambda}{4\epsilon} |1-v|^2 \f$ ). Dimension theory tells
    * that it is in 1/length unit, the lower the longer is the set of
    * discontinuities.
    */
    void setLambda( Scalar _lambda );

    /**
    * Sets the parameter \f$ \epsilon \f$ of AT functional. Should be
    * set \b after \ref setAlpha and \ref setLambda.
    *
    * @param _epsilon the \f$ \epsilon \f$ parameter in AT functional (
    * in terms \f$ \int \lambda \epsilon | v \grad u |^2 + \int
    * \frac{\lambda}{4\epsilon} |1-v|^2 \f$ ). Dimension theory tells
    * that it is in length unit, the lower the thinner is the set of
    * discontinuities.
    */
    void setEpsilon( Scalar _epsilon );

    /// Computes the SNR of u wrt ideal input (should have been given @see addInput).
    Scalar computeSNR() const;
    
    /// @return the (global) alpha parameter.
    Scalar getAlpha() const { return alpha; }

    /// @return the lambda parameter.
    Scalar getLambda() const { return lambda; }

    /// @return the epsilon parameter.
    Scalar getEpsilon() const { return epsilon; }

    /// @param i an integer (between 0 and the number of input forms).
    /// @return the \a i-th input \a g 0-form.
    const PrimalForm0& getG( int i ) const { return g0.at( i ); }

    /// @param i an integer (between 0 and the number of input forms).
    /// @return the \a i-th \a u 0-form.
    const PrimalForm0& getU( int i ) const { return u0.at( i ); }

    /// @return the \a v 1-form.
    const PrimalForm1& getV() const { return v1; }

    /// @return the size of a 0-form vector
    unsigned int size0() const { return alpha_Id0.myContainer.columns(); }

    /// @return the size of a 1-form vector
    unsigned int size1() const { return v1.myContainer.rows(); }

    // ----------------------- Solver --------------------------------------
  public:
    
    /// Computes a solution to function(s) \a u given the input \a g and current \a v.
    /// @return 'true' iff the solver worked.
    bool solveU();

    /// Computes a solution to function \a v given the input \a g and current \a u.
    /// @return 'true' iff the solver worked.
    bool solveV();

    /// Computes the variation of \a v after a call to \ref
    /// solveV. @see delta_v_l1, delta_v_l2, delta_v_loo.
    ///
    /// @return the max of all variations (i.e. delta_v_loo).
    Scalar computeVariation();

    /// Checks that form \a v is between 0 and 1 and forces \a v to be in-between.
    /// @return the max of the deviation wrt 0 and 1.
    Scalar checkV();

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;


    // ------------------------- Public Datas ------------------------------
  public:

    /// The L1-norm of variation of v.
    Scalar delta_v_l1;
    /// The L2-norm of variation of v.
    Scalar delta_v_l2;
    /// The Linfinity-norm of variation of v.
    Scalar delta_v_loo;

    /// edge laplacien
    PrimalIdentity1       L1;

    // ------------------------- Protected Datas ------------------------------
  protected:

    /// The g 0-forms
    std::vector< PrimalForm0 > g0;
    /// The ideal input 0-forms (for snr computation).
    std::vector< PrimalForm0 > i0;
    /// The u 0-forms
    std::vector< PrimalForm0 > u0;
    /// The v 1-form
    PrimalForm1 v1;
    /// The v 1-form at the previous iteration.
    PrimalForm1 former_v1;
    /// Smoothness parameter alpha of AT (in 1/area unit)
    double alpha;
    /// Amount of discontinuity parameter lambda (in 1/length unit).
    double lambda;
    /// Thickness of discontinuity set (in length unit).
    double epsilon;

    /// The solver for every 0-form u[i]
    SolverU solver_u;
    /// The solver for 1-form v
    SolverV solver_v;

    /// alpha Id0
    PrimalIdentity0 alpha_Id0;
    /// alpha g0
    std::vector< PrimalForm0 > alpha_g0;
    /// lambda * edge laplacien
    PrimalIdentity1       l_L1;
    /// lambda 1/4 1
    PrimalForm1           l_1_over_4;
    /// epsilon * lambda * edge laplacien + (lambda / (4*epsilon)) * Id1
    PrimalIdentity1       left_V1;
    /// lambda 1/(4*epsilon) 1
    PrimalForm1           l_1_over_4e;

    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Hidden services ------------------------------
  protected:

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ATu0v1


  /**
   * Overloads 'operator<<' for displaying objects of class 'ATu0v1'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ATu0v1' to write.
   * @return the output stream after the writing.
   */
  template <typename TKSpace, typename TLinearAlgebra>
  std::ostream&
  operator<< ( std::ostream & out, const ATu0v1<TKSpace, TLinearAlgebra> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "DGtal/dec/ATu0v1.ih"
#include "ATu0v1.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ATu0v1_h

#undef ATu0v1_RECURSES
#endif // else defined(ATu0v1_RECURSES)
