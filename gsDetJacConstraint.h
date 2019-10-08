/** @file gsDetJacConstraint.h

@brief  Class to handle constraints on the injectivity of a gsMultiPatch Parametrization
        It computes spline coefficients for the determinant of the Jacobian
        which then can be bounded to be always positive

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSDETJACCONSTRAINT_H
#define GSDETJACCONSTRAINT_H
#include "gsConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

class gsDetJacConstraint: public gsConstraint{
public:

    // Constructs from gsMultiPatch, runs the setup method.
    gsDetJacConstraint(gsMultiPatch<>* mpin, bool useTPSolver = false);

    // Get constraints (spline coefficients of detJ) by projection
    // FIXIT: Change the sign if the coefs is negative?
    gsVector<> evalCon(); // Implementation is in here..

    // Get derivatives
    // FIXIT: make dimension independent
    gsSparseMatrix<> getDerivRhsFromPatch(index_t patch);
    gsMatrix<> getJacobianFromPatch(index_t patch);
    gsIpOptSparseMatrix getJacobian();

    // Accessors
    const gsMultiBasis<> & detJacBasis() const { return m_detJacBasis; }
    const bool & isSolverSetup(index_t i) const { return m_areSolversSetup[i]; }

    // set the tolerance
    void setEps(real_t eps){m_eps = eps;}

    // Get the upper and lower bounds of constraint
    gsVector<> getUpperBounds();
    gsVector<> getLowerBounds();

    // Accessors
    index_t numConstraints(){ return m_size; };
    real_t eps(){ return m_eps; };

    // Get det J as a multipatch
    gsMultiPatch<> getDetJ();

    // Plot det J in paraview with filename "name"
    void plotDetJ(std::string name);

    // Plot active constraints
    // tol1 sets the tolerance for marking active constraints
    void plotActiveConstraints(std::vector<bool> & elMarked, std::string name, real_t tol1 = 0);

    // Get the mass matrix for patch i
    gsSparseMatrix<> getMassMatrix(index_t i);

    // Mark elements where constraints are active
    // tol1 sets the tolerance for marking active constraints
    // tol2 set the tolerance for which elements to refine, should be in [0,1]
    //      (if it is zero it refines support of active coefficients)
    void markElements(std::vector<bool> & elMarked, real_t tol1 = 0, real_t tol2 = 0);

    // Method to setup m_detJacBasis and m_space_mapper
    // It is called in the constructor and should be called whenever the space of m_mp change,
    //      (eg. when geometry is refined locally or globally)
    void setup();

    // Method to solve mass matrix system for arbitrary rhs
    gsVector<> massMatrixSolve(gsVector<> rhs) const;

    // Method to calculate the hessian of D_(ii), the right hand side for projection step
    gsSparseMatrix<> hessD(index_t patch, index_t i) const;

    index_t sizeOfBasis(index_t p){
        return m_detJacBasis.size(p);
    }

public:
    gsMultiBasis<> m_detJacBasis;

    real_t m_eps = 0;

    // use tensor product structure when solving system
    bool m_useTPSolver;
    std::vector<gsLinearOperator<>::Ptr> m_solversTensor;
    gsVector<gsSparseSolver<>::LU> m_solversMassMatrix;
    gsVector<bool> m_areSolversSetup;

    index_t m_size;

};




#endif //GSDETJACCONSTRAINT_H
