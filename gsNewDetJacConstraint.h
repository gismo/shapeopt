/** @file gsNewDetJacConstraint.h

@brief  Class to handle constraints on the injectivity of a gsMultiPatch Parametrization
        FIXIT: write explanation

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSNEWDETJACCONSTRAINT_H
#define GSNEWDETJACCONSTRAINT_H
#include "gsConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

class gsNewDetJacConstraint: public gsConstraint{
public:

    // Constructs from gsMultiPatch, runs the setup method.
    gsNewDetJacConstraint(memory::shared_ptr<gsMultiPatch<>> mpin);

    // Get constraints
    gsVector<> evalCon(); // Implementation is in here..

    // Get derivatives
    gsIpOptSparseMatrix getJacobian();

    // Set m_space_mapper, FIXIT clean up
    void setSpaceMapper();

    // set the tolerance
    void setEps(real_t eps){m_eps = eps;}

    // Get the upper and lower bounds of constraint
    gsVector<> getUpperBounds();
    gsVector<> getLowerBounds();

    // Accessors
    real_t eps(){ return m_eps; };

    // Get det J as a multipatch
    gsMultiPatch<> getDetJ();

    // Plot det J in paraview with filename "name"
    void plotDetJ(std::string name);

    // Plot active constraints
    // tol1 sets the tolerance for marking active constraints
    // FIXIT : implement this one
    void plotActiveConstraints(std::vector<bool> & elMarked, std::string name, real_t tol1 = 0);

    // Mark elements where constraints are active
    // tol1 sets the tolerance for marking active constraints
    // tol2 set the tolerance for which elements to refine, should be in [0,1]
    //      (if it is zero it refines support of active coefficients)
    // FIXIT: Implement this one
    void markElements(std::vector<bool> & elMarked, real_t tol1 = 0, real_t tol2 = 0)
    {
        GISMO_NO_IMPLEMENTATION;
    };

    // Method to setup method
    // It is called in the constructor and should be called whenever the space of m_mp change,
    //      (eg. when geometry is refined locally or globally)
    void setup();


public:
    real_t m_eps = 0;

    // Mappers for bases of partial derivatives
    std::vector< gsDofMapper > m_mappers;

    // Datastructure for storing "Overlapping" indicies
    gsMatrix<> m_overlapping;

};




#endif //GSNEWDETJACCONSTRAINT_H
