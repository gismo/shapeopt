#include <gismo.h>
#include "gsSpringMethod.h"
using namespace gismo;

// Debug
gsSpringMethod::gsSpringMethod(gsMultiPatch<> *mpin): gsAffineParamMethod(mpin),
    m_As(m_mp->targetDim()), m_bs(m_mp->targetDim()), m_solvers(m_mp->targetDim()){

    setupSolvers();
}

// Implement
gsSpringMethod::gsSpringMethod(gsMultiPatch<>* mpin,std::vector< gsDofMapper > mappers):
    gsAffineParamMethod(mpin,mappers), m_As(m_mp->targetDim()), m_bs(m_mp->targetDim()),
    m_solvers(m_mp->targetDim())
{
    setupSolvers();
}

void gsSpringMethod::setupSolvers()
{
    for (index_t d = 0; d < m_mp->targetDim(); d++){
        setupSystem(m_mappers[d],m_As[d],m_bs[d],d);
        m_solvers[d].compute(m_As[d]);
    }
}

// Implement,
// FIXIT make dimension independent
void gsSpringMethod::setupSystem(gsDofMapper &mapper, gsSparseMatrix<> &A, gsVector<> &b, index_t coord){
    index_t n = mapper.freeSize();
    // Prepare system and rhs
    gsSparseMatrix<> tmp(n,n);
    A.swap(tmp);
    A.reserve( gsVector<int>::Constant(A.cols(), 2*d) );
    b.setZero(n);

    // For each patch
    for (index_t p = 0; p < m_mp->nBoxes(); p++){
        gsVector<index_t,d> stride; // Tensor strides
        gsVector<index_t,d> size; // Tensor strides

        // Prepare stride and size
        static_cast<gsTensorBSplineBasis<d>&> (m_mp->patch(p).basis()).stride_cwise(stride);
        static_cast<gsTensorBSplineBasis<d>&> (m_mp->patch(p).basis()).size_cwise(size);

        const index_t dd = 2*d;

        for (int i = 0; i < m_mp->patch(p).coefsSize(); i++)
        {
            index_t ii = mapper.index(i,p);

            if ( mapper.is_coupled_index(ii)){
                if(mapper.is_free_index(ii)){
                    // gsInfo << "and is free\n";
                } else {
                    // gsInfo << "and is nonfree\n";
                }
            }

            if( !mapper.is_free_index(ii))
                continue;

            // if we are on the boundary of the gsMultiPatch we only have 3 neighbors!!
            // FIXIT: maybe use mappers instead of is_boundary
            if (is_boundary(i,p)){
                A(ii,ii) = dd - 1;
            } else {
                A(ii,ii) = dd;
            }

            for ( unsigned k = 0; k<d; k++ ) // for all neighbors
            {
                for ( int s = -1; s<2; s+=2 ) // +/- (up or down)
                {
                    if ( s==-1 && stride[k] == 1 && i % size[k] == 0 ) // If we are on left boundary we cannot decrease
                    {
                        continue;
                    }
                    if ( s==+1 && stride[k] == 1 && i % size[k] == size[k]-1 ) // If we are on right boundary we cannot increase
                    {
                        continue;
                    }

                    const unsigned j = i + s * stride[k];

                    if (j < 0 || j >= m_mp->patch(p).coefsSize()) // Out of bounds?
                    {
                        continue;
                    }

                    index_t jj = mapper.index(j,p);// get j-index in the matrix

                    if ( mapper.is_free_index(jj) )
                    {
                        A(ii,jj) = -1;
                        // gsInfo << ii << ", " << jj << " in A is altered...\n";
                    }
                    else // boundary node
                    {
                        // for 2D: add only half on the interfaces where this is traversed again later..
                        b[ii] += (mapper.is_coupled_index(ii) ? 0.5 : 1.0) * m_mp->patch(p).coef(j,coord);

                        // FIXIT think about whether this strategy is correct or not
                        // if (mapper.is_coupled_index(ii)) { // If we are on interface
                        //     // gsInfo << ii << " is coupled to ";
                        //
                        //     // Count no. of times the neighbor is coupled, this is the number of time we will visit it
                        //     std::vector< std::pair<index_t, index_t > > vec;
                        //     mapper.preImage(jj,vec);
                        //     real_t n_coupled = vec.size();
                        //
                        //     // Only add 1/n_coupled of the coef since it will be visited n_coupled times.
                        //     b[ii] += 1.0/n_coupled * m_mp->patch(p).coef(j,coord);
                        //
                        //     // gsInfo << vec.size() << "\n";
                        // } else { // If we are not on an interface, we can just add the neighbor to the rhs.
                        //     b[ii] +=  m_mp->patch(p).coef(j,coord);
                        // }
                    }
                }
            }
        }
    }
}

// Debug, maybe rethink?
gsVector<> gsSpringMethod::getUpdate(gsVector<> x)
{
    gsVector<> old_tagged = getTagged();

    updateTagged(x);
    gsVector<> out(n_free);
    for (index_t d = 0; d < m_mp->targetDim(); d++){
        gsSparseMatrix<> tmp;
        setupSystem(m_mappers[d], tmp, m_bs[d],d);
        out.segment(m_shift_free[d],m_mappers[d].freeSize()) = m_solvers[d].solve(m_bs[d]);
    }

    updateTagged(old_tagged);

    return out;
}

// FIXIT : really expensive method, preferably use the mapper instead
bool gsSpringMethod::is_boundary(index_t i, index_t p) const
{
    for (index_t k = 0; k < m_mp->nBoundary(); k ++){
        patchSide ps = m_mp->boundaries()[k];
        if (ps.patch != p) continue;

        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
        for (index_t j = 0; j < boundaryDofs.size(); j++){
            if (boundaryDofs[j] == i) return true;
        }
    }
    return false;
}
