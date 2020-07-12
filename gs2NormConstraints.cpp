#include <gismo.h>
#include "gsConstraint.h"
#include "gs2NormConstraints.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;


gs2NormConstraints::gs2NormConstraints(gsMultiPatch<>::Ptr mp, std::vector< gsDofMapper > mappers, real_t a, real_t b): 
    gsConstraint(mp),
    m_mappers(mappers),
    m_a(a), 
    m_b(b)
{

    //GISMO_ASSERT(m_mappers[0].taggedSize() == m_mappers[1].taggedSize(), "Error, the same cps need to be tagged for all directions");
    //GISMO_ASSERT(m_mappers[0].taggedSize() == m_mappers[2].taggedSize(), "Error, the same cps need to be tagged for all directions");

    m_dim = m_mp->targetDim();

    n_free = 0;
    m_free_shift.setZero(m_dim);;
    for (index_t d = 0; d < m_dim; d++)
    {
        n_free += m_mappers[d].freeSize();
        if ( d > 0)
            m_free_shift[d] = m_free_shift[d-1] + m_mappers[d-1].freeSize();
    }

    gsDebugVar(m_free_shift);

    // 
    n_cps = m_mappers[0].size();
    gsDebugVar(n_cps);
    gsDebugVar(n_free);

    // Initialize constraintMap
    constraintMap.reserve(n_cps);
    for( index_t ii = 0; ii < n_cps; ii++)
        constraintMap[ii] = -1;

    // Calculate n_constraints and setup constraint map
    index_t k = 0; // constraint counter

    // For each cps component
    for( index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        for (index_t i = 0; i < m_mp->patch(p).coefsSize(); i++)
        {
            // Is there a component of (i,p) that are free 
            //     and is there a component of (i,p) that are tagged?
            //gsInfo << "(" << i << ", " << p << ") : is_free = " << isACompFree(i,p) << " , is_tagged = " << isACompTagged(i,p) << "\n";
            if ( isACompFree(i,p) && isACompTagged(i,p) ) 
            {
                index_t ii = index(i,p);
                if ( constraintMap[ii] == -1) // If we did not flag ii
                {
                    constraintMap[ii] = k;
                    k++;
                }
            }
        }
    }
    n_constraints = k;
    gsDebugVar(n_constraints);

}

gsVector<> gs2NormConstraints::evalCon() 
{
    gsVector<> out(n_constraints);

    for (index_t ii = 0; ii < n_cps; ii++)
    {
        index_t k = constraintMap[ii];

        if (k != -1) // If a constraint is attached to ii
        {
            // get local index (p,i);
            std::pair< index_t, index_t > local = lindex(ii); // WHAT SHOULD d BE?
            index_t p = local.first;
            index_t i = local.second;

            out[k] = m_mp->patch(p).coef(i).squaredNorm();
            //gsDebugVar(out[k]);
        }
    }

    return out;

}

gsIpOptSparseMatrix gs2NormConstraints::getJacobian() 
{

    gsMatrix<> jac(n_constraints,n_free);
    jac.setZero(n_constraints,n_free);

    for (index_t ii = 0; ii < n_cps; ii++)
    {
        index_t k = constraintMap[ii];

        if (k != -1) // If a constraint is attached to ii
        {
            // get local index (p,i) using m_mappers[0] index
            std::pair< index_t, index_t > local = lindex(ii); // We use m_mappers[0] to map 
            index_t p = local.first;
            index_t i = local.second;

            for (index_t d = 0; d < m_dim; d++) // For each component
            {
                // Get the m_mappers[d] index
                if (!isFree(i,p,d)) continue; // Skip if i,p,d is not free

                index_t jj = index(i, p, d); // Get global index of (i,p,d)
                jac(k,jj) = 2*m_mp->patch(p).coef(i,d);

            }
        }
    }

    // FIXIT why are the norm sometimes 0?
    //gsDebugVar(jac.norm());

    // We have to use a dense matrix otherwise I got an error.! Perhaps fix this? FIXIT
    gsIpOptSparseMatrix out(jac, -1); 

    return out;
}

gsIpOptSparseMatrix gs2NormConstraints::getJacobian(gsVector<> free) 
{

    gsMatrix<> jac(n_constraints,n_free);
    jac.setZero(n_constraints,n_free);

    for (index_t ii = 0; ii < n_cps; ii++)
    {
        index_t k = constraintMap[ii];

        if (k != -1) // If a constraint is attached to ii
        {
            // get local index (p,i) using m_mappers[0] index
            std::pair< index_t, index_t > local = lindex(ii); // We use m_mappers[0] to map 
            index_t p = local.first;
            index_t i = local.second;

            for (index_t d = 0; d < m_dim; d++) // For each component
            {
                // Get the m_mappers[d] index
                if (!isFree(i,p,d)) continue; // Skip if i,p,d is not free

                index_t jj = index(i, p, d); // Get global index of (i,p,d)
                jac(k,jj) = 2*m_mp->patch(p).coef(i,d);

                gsInfo << "const " << k << ": " << i << ", " << p << " attach to " << jj << "\n"; 

                gsDebugVar(m_mp->patch(p).coef(i,d) - free[jj]);
            }
        }
    }

    // FIXIT why are the norm sometimes 0?
    //gsDebugVar(jac.norm());

    // We have to use a dense matrix otherwise I got an error.! Perhaps fix this? FIXIT
    gsIpOptSparseMatrix out(jac, -1); 

    return out;
}

gsVector<> gs2NormConstraints::getUpperBounds() 
{
    gsVector<> out(n_constraints);

    for (index_t i = 0; i < n_constraints; i++)
    {
        out[i] = m_b*m_b;
    }

    return out;

}

gsVector<> gs2NormConstraints::getLowerBounds() 
{
    gsVector<> out(n_constraints);

    for (index_t i = 0; i < n_constraints; i++)
    {
        out[i] = m_a*m_a;
    }

    return out;

}

bool gs2NormConstraints::isACompTagged( index_t i , index_t p)
{
    for (index_t d = 0; d < m_dim; d++)
    {
        if (m_mappers[d].is_tagged(i,p)) // If (i,p,d) is tagged
            return true;
    }
    return false; // We did not find a tagged one

}

bool gs2NormConstraints::isFree( index_t i, index_t p, index_t d)
{
    return m_mappers[d].is_free(i,p); // If (i,p,d) is free
}

bool gs2NormConstraints::isACompFree( index_t i, index_t p)
{
    for (index_t d = 0; d < m_dim; d++)
    {
        if (m_mappers[d].is_free(i,p)) // If (i,p,d) is free
            return true;
    }
    return false; // We did not find a tagged one

}

index_t gs2NormConstraints::index( index_t i, index_t p, index_t d )
{
    return m_free_shift[d] + m_mappers[d].index(i,p);
}

std::pair< index_t, index_t > gs2NormConstraints::lindex( index_t ii, index_t d )
{
    std::vector< std::pair<index_t,index_t> > result;
    m_mappers[d].preImage(ii, result);

    return result[0];
}
