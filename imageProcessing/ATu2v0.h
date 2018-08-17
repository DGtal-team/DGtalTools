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
 * @file ATu2v0.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Marion Foare (\c marion.foare@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2016/10/12
 *
 * Header file for module ATu2v0.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ATu2v0_RECURSES)
#error Recursive header files inclusion detected in ATu2v0.h
#else // defined(ATu2v0_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ATu2v0_RECURSES

#if !defined ATu2v0_h
/** Prevents repeated inclusion of headers. */
#define ATu2v0_h

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
  // template class ATu2v0
  /**
  * Description of template class 'ATu2v0' <p> \brief Aim: This class
  * solves Ambrosio-Tortorelli functional in a plane for \a u a
  * (vector of) 2-form(s) and \a v a 0-form. \a u is a regularized
  * approximation of an input image data \a g (grey-level or color
  * image), while \a v represents the set of discontinuities of \a u.
  *
  * @tparam TKSpace any model of CCellularGridSpaceND, e.g KhalimskySpaceND
  * @tparam TLinearAlgebra any back-end for performing linear algebra, default is EigenLinearAlgebraBackend.
  *
  */
  template < typename TKSpace,
             typename TLinearAlgebra = EigenLinearAlgebraBackend >
  struct ATu2v0 : public DECImage2D<TKSpace, TLinearAlgebra>
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
    using typename Base::PrimalAntiderivative1;
    using typename Base::PrimalAntiderivative2;
    using typename Base::DualAntiderivative1;
    using typename Base::DualAntiderivative2;
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
    typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 2, PRIMAL, 2, PRIMAL> 
                                                           SolverU;
    typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 0, PRIMAL, 0, PRIMAL> 
                                                           SolverV;

    BOOST_STATIC_ASSERT(( KSpace::dimension == 2 ));

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~ATu2v0() = default;

    /**
    * Default constructor. The object needs to be initialized with \ref init.
    * @param _verbose specifies the verbose level (0: silent, 1: more info ... ). 
    */
    ATu2v0( int _verbose = 1 );
    
    /**
    * Constructor from Khalimsky space, which specifies the domain of calculus.
    */
    void init( Clone<KSpace> K );

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    ATu2v0 ( const ATu2v0 & other ) = delete;

    /**
     * Move constructor.
     * @param other the object to move.
     */
    ATu2v0 ( ATu2v0 && other ) = delete;

    /**
     * Copy assignment operator.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    ATu2v0 & operator= ( const ATu2v0 & other ) = delete;

    /**
     * Move assignment operator.
     * @param other the object to move.
     * @return a reference on 'this'.
     */
    ATu2v0 & operator= ( ATu2v0 && other ) = delete;

    /**
    * Adds an input 2-form by filtering an \a image values.
    *
    * @param image any image such that the domain of this space is
    * included in the domain of the image.
    *
    * @param f any functor associated a scalar to an image value.
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
    template <typename Image, typename Function>
    void addInput( const Image& image,
                   const Function& f,
                   bool perfect_data = false );

    /// Use metric average to smooth L1-metric effects.
    void setMetricAverage( bool average );
    
    /// Sets approximation \a u to be equal to the input. Used for
    /// initializating \a u. Should be called once all \ref addInput
    /// have been called.
    void setUFromInput();

    /// Sets approximation \a u to be equal to the input. Used for
    /// initializating \a u. Should be called once all \ref addInput
    /// have been called. Note that it initializes \a u with random
    /// values wherever an inpainting mask was applied.
    /// @see setAlpha( Scalar _alpha, const PrimalForm2& m )
    void setUFromInputAndMask();

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
    void setAlpha( Scalar _alpha, const PrimalForm2& m );

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
    /// @return the \a i-th input \a g 2-form.
    const PrimalForm2& getG( int i ) const { return g2.at( i ); }

    /// @param i an integer (between 0 and the number of input forms).
    /// @return the \a i-th \a u 2-form.
    const PrimalForm2& getU( int i ) const { return u2.at( i ); }

    /// @return the \a v 0-form.
    const PrimalForm0& getV() const { return v0; }

    /// @return the size of a 0-form vector
    unsigned int size0() const { return v0.myContainer.rows(); }

    // /// @return the size of a 1-form vector
    // unsigned int size1() const { return v1.myContainer.rows(); }

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

    /// point_to_edge average operator
    PrimalDerivative0     M01;
    /// edge_to_face average operator
    PrimalDerivative1     M12;
    /// Antiderivative 2-form -> to 1-form
    PrimalAntiderivative2 primal_AD2;

    // ------------------------- Protected Datas ------------------------------
  protected:

    /// The g 2-forms
    std::vector< PrimalForm2 > g2;
    /// The ideal input 2-forms (for snr computation).
    std::vector< PrimalForm2 > i2;
    /// The u 2-forms
    std::vector< PrimalForm2 > u2;
    /// The v 0-form
    PrimalForm0 v0;
    /// The v 0-form at the previous iteration.
    PrimalForm0 former_v0;
    /// Smoothness parameter alpha of AT (in 1/area unit)
    double alpha;
    /// Amount of discontinuity parameter lambda (in 1/length unit).
    double lambda;
    /// Thickness of discontinuity set (in length unit).
    double epsilon;

    /// The solver for every 2-form u[i]
    SolverU solver_u;
    /// The solver for 0-form v
    SolverV solver_v;

    /// alpha Id2
    PrimalIdentity2 alpha_Id2;
    /// alpha g2
    std::vector< PrimalForm2 > alpha_g2;
    /// (lambda / (4*epsilon)) * Id0 + lambda epsilon D0^t D0
    PrimalIdentity0       left_V0;
    /// lambda 1/4 1
    PrimalForm0           l_1_over_4;
    /// lambda 1/(4*epsilon) 1
    PrimalForm0           l_1_over_4e;
    /// lambda 1/(4_epsilon) Id0
    PrimalIdentity0       l_1_over_4e_Id0;    

    /// When 'true', use metric average, otherwise use identity.
    bool metric_average;
    
    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Hidden services ------------------------------
  protected:

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ATu2v0


  /**
   * Overloads 'operator<<' for displaying objects of class 'ATu2v0'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ATu2v0' to write.
   * @return the output stream after the writing.
   */
  template <typename TKSpace, typename TLinearAlgebra>
  std::ostream&
  operator<< ( std::ostream & out, const ATu2v0<TKSpace, TLinearAlgebra> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "DGtal/dec/ATu2v0.ih"
#include "ATu2v0.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ATu2v0_h

#undef ATu2v0_RECURSES
#endif // else defined(ATu2v0_RECURSES)
