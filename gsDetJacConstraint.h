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
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

class gsDetJacConstraint{
public:

    // Constructs from gsMultiPatch
    gsDetJacConstraint(gsMultiPatch<>* mpin);

    // Get constraints (spline coefficients of detJ) by projection
    // FIXIT: Change the sign if the coefs is negative?
    void evalCon_into(gsAsVector<real_t> & result);
    gsVector<> evalCon(); // Implementation is in here..

    // Get derivatives
    // FIXIT: make dimension independent
    void getDerivRhsFromPatch(index_t patch, gsSparseMatrix<> &xJac, gsSparseMatrix<> &yJac);
    void getJacobianFromPatch(index_t patch, gsMatrix<> &xJac, gsMatrix<> &yJac);
    gsIpOptSparseMatrix getJacobian();
    void jacobCon_into(gsAsVector<real_t> & result);

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

    void plotDetJ(std::string name);
public:
    gsMultiPatch<>* m_mp;
    gsMultiBasis<> m_detJacBasis;

    real_t m_eps;

    gsVector<gsSparseSolver<>::LU> m_solversMassMatrix;
    gsVector<bool> m_areSolversSetup;

    index_t m_size;
    index_t n_controlpoints;
    index_t n_constraints;

};




#endif //GSDETJACCONSTRAINT_H
