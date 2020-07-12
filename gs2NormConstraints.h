/** @file gs2NormConstraint.h

@brief  Base class to handle constraints on the two norm of tagged controlpoints

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde
*/
#ifndef GS2NORMCONSTRAINTS_H
#define GS2NORMCONSTRAINTS_H
#include "gsConstraint.h"
using namespace gismo;

class gs2NormConstraints : public gsConstraint {
public:
    gs2NormConstraints(gsMultiPatch<>::Ptr mp, std::vector< gsDofMapper > mappers, real_t a, real_t b);

    gsVector<> evalCon();

    gsIpOptSparseMatrix getJacobian();
    gsIpOptSparseMatrix getJacobian(gsVector<> free);

    gsVector<> getUpperBounds();
    gsVector<> getLowerBounds();

    // Helper methods
    bool isACompTagged( index_t i, index_t p);    // Is local index tagged in some component
    bool isACompFree( index_t i, index_t p);      // Is local index free in some component
    bool isFree( index_t i, index_t p, index_t d);      // Is local index (i,p,d) free ?
    index_t index( index_t i, index_t p, index_t d = 0);       // Get global index from (i,p,d)
    std::pair< index_t, index_t > lindex( index_t ii, index_t d = 0);                     // Get local index (p,i) from global

public:

    typedef memory::shared_ptr< gs2NormConstraints > Ptr;
    typedef memory::unique_ptr< gs2NormConstraints > uPtr;

public:
    
    std::vector< gsDofMapper > m_mappers;

    std::vector< index_t > constraintMap; // Map from cps (i,p) to constraint. Initialized to -1.

    index_t n_cps;
    index_t n_free;

    gsVector< index_t > m_free_shift; // Shifts of cps for each component. Used on index(i,p,d).

    real_t m_a, m_b;

    index_t m_dim;

};




#endif //GS2NORMCONSTRAINTS_H
