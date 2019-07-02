/** @file gsIpOptSparseMatrix.h

@brief  A Implements a base class for parametrization method for multipatch geometry
It uses a list of gsDofMappers (one for each coordinate) which to denote wich
DoFs that are fixed and which are not fixed.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSPARAMMETHOD_H
#define GSPARAMMETHOD_H
using namespace gismo;

#include "gsIpOptSparseMatrix.h"

class gsParamMethod{
public:

    // Constructs from multipatch, should setup mappers and assume interfaces glued
    // and boundary fixed; FIXIT: implement this functionality, FIXIT: tag bnd?
    gsParamMethod(gsMultiPatch<>* mpin);

    // Constructs from list of mappers
    gsParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers);

    // Sets the mappers and (re-)computes n_free, n_tagged, n_cps, and shifts
    void setMappers(std::vector< gsDofMapper > mappers);

    // Methods to reset the parametrization method
    // Currently only used for the gsAffineOptParamMethod, but it might be needed in other classes as well.
    // The default behaviour is to do nothing.
    virtual void reset(){};

    virtual void updateAndReset() { GISMO_NO_IMPLEMENTATION; };

    // FIXIT: Think about which one of these functions is the best to have virtual...
    //        i.e. one could be implemented generic by the other one...
    // Method to update inner controlpoints. Should be implemented in inherited classes
    virtual void update(){ GISMO_NO_IMPLEMENTATION; };

    // Method to update inner controlpoints, given tagged Dofs (x).
    // Should be implemented in inherited classes
    virtual void update(gsVector<> x){ GISMO_NO_IMPLEMENTATION; };

    // Method to get the jacobian of the update (free Dofs) with respect to x (the tagged Dofs)
    // Should be implemented in inherited classes
    virtual gsMatrix<> jacobUpdate(gsVector<> x){ GISMO_NO_IMPLEMENTATION; };

    // Get a vector of the tagged control points.
    gsVector<> getTagged();

    // Update the tagged control points.
    void updateTagged(gsVector<> x) const;

    // Returns vector with free variables (with coordinates stacked, first x then y the possibly z)
    gsVector<> getFree() const;

    // Update the gsMultiPatch* m_mp, based on a vector of free variables
    void updateFree(gsVector<> des) const;

    // Update the gsMultiPatch* m_mp, based on a vector of free variables, and vector of tagged cps
    void updateFreeAndTagged(gsVector<> des, gsVector<> x) const;

    // Get the vector of all control points (eliminated and free)
    gsVector<> getControlPoints() const;

    // Get the vector of all control points (eliminated and free) from design variables
    gsVector<> getControlPoints(gsVector<> des) const;

    // Updates m_mp by vector of all controlpoints (e.g. obtained from getControlPoints() )
    // OBS: also updates eliminated DoFs...
    void updateControlPoints(gsVector<> cps);

    // Get the vector of all control points in flat layout (one patch at a time)
    gsVector<> getFlat() const;

    // Update from vector of all control points in flat layout (one patch at a time)
    void updateFlat(gsVector<> flat) const;

    // Setup the default mappers, ie. fix and tag the boundary and free the rest (including interfaces)
    void setupDefaultMappers();

    // Method overloaded by gsAffineParamMethod
    virtual void computeMap(){ GISMO_NO_IMPLEMENTATION; };

    // Accessors
    std::vector< gsDofMapper > mappers() { return m_mappers; };
    gsMultiPatch<>* mp() { return m_mp; };

    // Maps a vector from mapper_in indexing, to m_mappers.
    // E.g. used to map gradients with respect to all cps to be respect to free
    // DoFs...
    // FIXIT: maybe find a more elegant way to handle all of this..
    //      Maybe always assemble wrt. specific mapper, and hold a permutation
    //      in a matrix, that we only need to multiply to map...
    gsMatrix<> mapMatrix(gsDofMapper mapper_in, gsMatrix<> mat) const;
    gsMatrix<> mapMatrixToTagged(gsDofMapper mapper_in, gsMatrix<> mat) const;

    // Maps a gsIpOptSparseMatrix from mapper_in indexing, to m_mappers.
    // E.g. used to map gradients with respect to all cps to be respect to free
    // DoFs...
    gsIpOptSparseMatrix mapMatrix(gsDofMapper mapper_in, gsIpOptSparseMatrix M) const;

public:
    mutable gsMultiPatch<>* m_mp;

    std::vector< gsDofMapper > m_mappers; // Mapper for each coordinate

    gsVector<> m_shift_free; // Contains the shifts used when mapping free DoFs
    gsVector<> m_shift_all; // Contains the shifts used when mapping all DoFs

    index_t n_free; // sum over coordinates
    index_t n_cps; // sum over coordinates
    index_t n_flat; // sum over coordinates
    index_t n_tagged; // sum over coordinates
    // index_t n_controlpoints;
    index_t fixedPatch = 3; // FIXIT: Should be moved to gsOptAntenna.

};

# endif //GSPARAMMETHOD_H
