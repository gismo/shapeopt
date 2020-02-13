/** @file gsConstraint.h

@brief  Base class to handle constraints, which can be used for optimization with IpOpt
        in a seperate class

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSCONSTRAINT_H
#define GSCONSTRAINT_H
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

class gsConstraint{
public:

    gsConstraint(){}; // Empty constructor

    // Constructs from gsMultiPatch
    gsConstraint(memory::shared_ptr<gsMultiPatch<>> mpin);

    // Evaluate into vector, currently calls evalCon and copy, but can be overloaded to improve efficiency.
    virtual void evalCon_into(gsAsVector<real_t> & result);

    // Evaluate constraint, should be overloaded in derived class
    virtual gsVector<> evalCon(){ GISMO_NO_IMPLEMENTATION; };

    // Evaluate derivatives, should be overloaded in derived class
    virtual gsIpOptSparseMatrix getJacobian(){ GISMO_NO_IMPLEMENTATION; };

    // Access mapper for space used for assembling derivatives
    // FIXIT, find more elegant way, when derivativs do not depend on derivatives being assembled
    const gsDofMapper & space_mapper() const { return m_space_mapper; }

    // Get the upper and lower bounds of constraint
    virtual gsVector<> getUpperBounds(){ GISMO_NO_IMPLEMENTATION; };
    virtual gsVector<> getLowerBounds(){ GISMO_NO_IMPLEMENTATION; };

    // Accessors
    index_t numConstraints(){ return n_constraints; };

    // Get the sign of detJ of a patch
    index_t getSignOfPatch(index_t patch);

    // Method to call after reseting the space of mpin, for example after refinement
    virtual void setup(){};

    // Method to mark elements for refinement
    virtual void markElements(std::vector<bool> & elMarked, real_t tol1 = 0, real_t tol2 = 0)
    {
        GISMO_NO_IMPLEMENTATION;
    };

    // Method to solve mass matrix system for arbitrary rhs
    virtual gsVector<> massMatrixSolve(gsVector<> rhs) const
    {
        GISMO_NO_IMPLEMENTATION;
    };


    // Method to calculate the hessian of D_(ii), the right hand side for projection step
    virtual gsSparseMatrix<> hessD(index_t patch, index_t i) const
    {
        GISMO_NO_IMPLEMENTATION;
    };


    virtual index_t sizeOfBasis(index_t p)
    {
        GISMO_NO_IMPLEMENTATION;
    };

    // Method to for example evaluate constraint that is getting aggregated
    virtual gsVector<> evalSubCon()
    {
        GISMO_NO_IMPLEMENTATION;
    };

    virtual real_t refineUntilPositive(index_t maxRefSteps = 10, real_t tol = 0)
    {
        GISMO_NO_IMPLEMENTATION;
    };

    void setEps(real_t eps) { m_eps = eps; };

public:
    memory::shared_ptr<gsMultiPatch<>> m_mp;

    // Holds the mapper used to compute the jacobian of the constraint
    // for extracting DoFs later
    // FIXIT: find prettier solution, since derivative might not be calculated using gsExprAssembler
    gsDofMapper m_space_mapper;

    index_t n_controlpoints;
    index_t n_constraints;

    real_t m_eps = 0;

};




#endif //GSCONSTRAINT_H
