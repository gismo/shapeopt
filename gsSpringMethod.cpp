#include <gismo.h>
#include "springMethod.h"
using namespace gismo;

springMethod::springMethod(gsMultiPatch<> *mpin): paraWithMapper(mpin), dJC(mpin), A(mp->targetDim()), b(mp->targetDim()), solvers(mp->targetDim()){
    setupMapper();
    ndesign = m_mappers[0].freeSize() + m_mappers[1].freeSize();
    setupSystem(m_mappers[0],A[0],b[0],0);
    setupSystem(m_mappers[1],A[1],b[1],1);
    // gsInfo << static_cast<gsMatrix<>>(Ax) << "\n";
    // gsInfo << bx << "\n";
    solvers[0].compute(A[0]);
    solvers[1].compute(A[1]);

    refCps = getControlPoints();
}

springMethod::springMethod(paraOptProblem* problem): paraWithMapper(problem->mp), dJC(problem->mp){
    setupMapper();
    ndesign = m_mappers[0].freeSize() + m_mappers[1].freeSize();
    setupSystem(m_mappers[0],A[0],b[0],0);
    setupSystem(m_mappers[1],A[1],b[1],1);
    // gsInfo << static_cast<gsMatrix<>>(Ax) << "\n";
    // gsInfo << bx << "\n";
    solvers[0].compute(A[0]);
    solvers[1].compute(A[1]);

    refCps = getControlPoints();
}

void springMethod::setupSystem(gsDofMapper &mapper, gsSparseMatrix<> &A, gsVector<> &b, index_t coord){
    index_t n = mapper.freeSize();
    // Prepare system and rhs
    gsSparseMatrix<> tmp(n,n);
    A.swap(tmp);
    A.reserve( gsVector<int>::Constant(A.cols(), 2*d) );
    b.setZero(n);

    // For each patch
    for (index_t p = 0; p < mp->nBoxes(); p++){
        gsVector<index_t,d> stride; // Tensor strides
        gsVector<index_t,d> size; // Tensor strides

        // Prepare stride and size
        static_cast<gsTensorBSplineBasis<d>&> (mp->patch(p).basis()).stride_cwise(stride);
        static_cast<gsTensorBSplineBasis<d>&> (mp->patch(p).basis()).size_cwise(size);

        const index_t dd = 2*d;

        for (int i = 0; i < mp->patch(p).coefsSize(); i++)
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
            if (is_boundary(i,p)){
                A(ii - mapper.firstIndex(),ii- mapper.firstIndex()) = dd - 1;
            } else {
                A(ii- mapper.firstIndex(),ii- mapper.firstIndex()) = dd;
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

                    if (j < 0 || j >= mp->patch(p).coefsSize()) // Out of bounds?
                    {
                        continue;
                    }

                    index_t jj = mapper.index(j,p);// get j-index in the matrix

                    if ( mapper.is_free_index(jj) )
                    {
                        jj -= mapper.firstIndex(); // Substract firstIndex to get rid of any possible shift
                        A(ii - mapper.firstIndex(),jj) = -1;
                        // gsInfo << ii << ", " << jj << " in A is altered...\n";
                    }
                    else // boundary node
                    {
                        // for 2D: add only half on the interfaces where this is traversed again later..
                        b[ii - mapper.firstIndex()] += (mapper.is_coupled_index(ii) ? 0.5 : 1.0) * mp->patch(p).coef(j,coord);

                        // FIXIT: only once for pair ii,j ???
                    }
                }
            }
        }
    }
}

gsVector<> springMethod::solveSystems(gsVector<> deltaCps){
    // update control points FIXIT should be improved to a more flexible way..
    dJC.updateDesignVariables(deltaCps);  // Set controlpoints to deltaCps, temporarily that is

    // DO SOMETHING TO MATCH INTERFACES !!!
    // Save nonzero values
    gsVector<> globalVec;
    globalVec.setZero(m_mappers[0].size() + m_mappers[1].size());
    for(index_t p = 0; p < mp->nBoxes(); p++){
        for (index_t i = 0; i < mp->patch(p).coefsSize(); i ++){
            if (mp->patch(p).coef(i,0) != 0){
                index_t iix = m_mappers[0].index(i,p);
                globalVec[iix] = mp->patch(p).coef(i,0);
            }
            if (mp->patch(p).coef(i,1) != 0){
                index_t iiy = m_mappers[1].index(i,p);
                // We need to correct the index since it is shifted with freeSize of m_mappers[0]
                globalVec[iiy - m_mappers[1].firstIndex() + m_mappers[0].size()] = mp->patch(p).coef(i,1);
            }
        }
    }

    dJC.updateDesignVariables(refCps);  // Set controlpoints to refCps
    // Set all control points by adding global vec
    for(index_t p = 0; p < mp->nBoxes(); p++){
        for (index_t i = 0; i < mp->patch(p).coefsSize(); i ++){
            index_t iix = m_mappers[0].index(i,p);
            mp->patch(p).coef(i,0) += globalVec[iix];
            index_t iiy = m_mappers[1].index(i,p);
                // We need to correct the index since it is shifted with freeSize of m_mappers[0]
            mp->patch(p).coef(i,1) += globalVec[iiy - m_mappers[1].firstIndex() + m_mappers[0].size()];
        }
    }

    // Get the new rhs vectors
    gsSparseMatrix<> tmp;
    setupSystem(m_mappers[0],tmp,b[0],0);
    setupSystem(m_mappers[1],tmp,b[1],1);

    // Get new design variables
    gsVector<> des(ndesign);
    gsVector<> des_x = solvers[0].solve(b[0]);
    gsVector<> des_y = solvers[1].solve(b[1]);
    des << des_x,des_y;

    updateDesignVariables(des);
    return des;
}

// Right now it returns a vector with all the control points
gsVector<> springMethod::solve(gsVector<> deltaCps){
    gsVector<> oldCps = dJC.getDesignVariables();
    gsVector<> out = getControlPoints(solveSystems(deltaCps)) - refCps;
    dJC.updateDesignVariables(oldCps);
    return out;
};

void springMethod::solve(){
    gsVector<> deltaCps = dJC.getDesignVariables();
    deltaCps.setZero();
    solveSystems(deltaCps);
};

void springMethod::solveAndUpdate(gsVector<> deltaCps){
    gsVector<> des = solveSystems(deltaCps);
    // updateDesignVariables(des);

};

bool springMethod::is_boundary(index_t i, index_t p) const {
    // FIXIT : quite expensive method
    for (index_t k = 0; k < mp->nBoundary(); k ++){
        patchSide ps = mp->boundaries()[k];
        if (ps.patch != p) continue;

        gsVector<unsigned> boundaryDofs = mp->basis(ps.patch).boundary(ps);
        for (index_t j = 0; j < boundaryDofs.size(); j++){
            if (boundaryDofs[j] == i) return true;
        }
    }
    return false;
}
