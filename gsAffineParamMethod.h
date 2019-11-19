/** @file gsAffineParamMethod.h

@brief  Implements an affine parametrization method. The main difference from
       gsParamMethod is that since we know that the parametrization method is
       affine, we precompute its value at 0 and its matrix, with respect to
       the tagged DoF's in the mappers, for efficient computations

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSAFFINEPARAMMETHOD_H
#define GSAFFINEPARAMMETHOD_H
using namespace gismo;

#include "gsParamMethod.h"

class gsAffineParamMethod: public gsParamMethod{
public:
    // Empty constructor
    gsAffineParamMethod(){ };

    // Constructs from multipatch, by eliminating boundary, gluing interfaces and tagging bnd
    gsAffineParamMethod(gsMultiPatch<>* mpin);

    // Constructs from mappers, one for each coordinate, should be finalized with
    // design variables for shape optimization tagged.
    gsAffineParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers);

    // Update controlpoints given by this method
    // OBS, as of now, computeMap has to be called first!
    // Returns false if the update failed
    bool update(gsVector<> x); // Ax + b

    // Update controlpoints, does not use A and b, but uses getUpdate
    // Returns false if the update failed
    bool update();

    // Derivatives of update with respect to tagged Dofs
    // OBS, as of now, computeMap has to be called first!
    gsMatrix<> jacobUpdate(gsVector<> x); // A

    // Method to get update as a vector, to be used to calculate A and b
    // Depends on the specific parametrization, ie. should be overloaded in
    // inherited class
    // Returns: vector of free control points
    virtual gsVector<> getUpdate(gsVector<> x) = 0;

    // Method to compute A and b
    void computeMap();


public:
    gsMatrix<> A; // Matrix that defines the affine map
    gsVector<> b; // Value at 0

};

#endif //GSAFFINEPARAMMETHOD_H
