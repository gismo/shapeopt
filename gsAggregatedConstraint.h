/** @file gsAggregatedConstraint.h

@brief  Class to aggregate constraints. The idea is to compose a constraint f:R_d - R_c with g:R_c -> R_{tilde{c}}
        with R_d the design space, R_c the space for constraints and R_{tilde{c}} the new constraint space.
        One reason to do this is to degrease the number of constraints, and therefore the size of the jacobian.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSAGGREGATEDCONSTRAINT_H
#define GSAGGREGATEDCONSTRAINT_H
#include "gsConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

class gsAggregatedConstraint: public gsConstraint{
public:

    // Constructs from gsConstraint and gsFunction, runs the setup method.
    // \a constraint will be the constraint
    // \a map is the function the you want to compose with
    gsAggregatedConstraint(gsMultiPatch<>* mp, gsConstraint* constraint, gsFunction<>* map);

    // Get constraints
    gsVector<> evalCon();

    // Get derivatives
    gsIpOptSparseMatrix getJacobian();

    // Set m_space_mapper, DO WE NEED THIS ONE?
    // void setSpaceMapper();

    // set the tolerance
    void setEps(real_t eps){m_eps = eps;}

    // Get the upper and lower bounds of constraint
    gsVector<> getUpperBounds();
    gsVector<> getLowerBounds();

    // Accessors
    real_t eps(){ return m_eps; };

    // Method to setup method
    // It is called in the constructor and should be called whenever the space of m_mp change,
    //      (eg. when geometry is refined locally or globally)
    void setup();

    gsVector<> evalSubCon();


public:
    real_t m_eps = 0;

    gsMultiPatch<>* m_mp;

    gsConstraint*   m_constraint;
    gsFunction<>*   m_map;

};




#endif //GSAGGREGATEDCONSTRAINT_H
