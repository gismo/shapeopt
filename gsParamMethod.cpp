#include <gismo.h>
#include <fstream>
#include "gsParamMethod.h"
using namespace gismo;

gsParamMethod::gsParamMethod(gsMultiPatch<>* mpin,std::vector< gsDofMapper > mappers):
m_mp(mpin),
m_isBoundaryFixed(mpin->targetDim(),mpin->nBoundary()),
m_isInterfaceFixed(mpin->targetDim(),mpin->nInterfaces()),
m_isBoundaryTagged(mpin->targetDim(),mpin->nBoundary()),
m_isInterfaceTagged(mpin->targetDim(),mpin->nInterfaces())
{
    setMappers(mappers);
}

gsParamMethod::gsParamMethod(gsMultiPatch<>* mpin):
m_mp(mpin),
m_mappers(m_mp->targetDim()),
m_isBoundaryFixed(mpin->targetDim(),mpin->nBoundary()),
m_isInterfaceFixed(mpin->targetDim(),mpin->nInterfaces()),
m_isBoundaryTagged(mpin->targetDim(),mpin->nBoundary()),
m_isInterfaceTagged(mpin->targetDim(),mpin->nInterfaces())
{
    // Setup default mappers
    setupDefaultMappers();
}

void gsParamMethod::setMappers(std::vector< gsDofMapper > mappers)
{
    m_mappers = mappers;

    m_shift_free.setZero(m_mp->targetDim());
    m_shift_all.setZero(m_mp->targetDim());

    n_free = 0;
    n_cps = 0;
    n_tagged = 0;
    n_flat = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_free += m_mappers[i].freeSize();
        n_cps += m_mappers[i].size();
        n_tagged += m_mappers[i].taggedSize();

        for(index_t p = 0; p < m_mp->nBoxes(); p++){
            n_flat += m_mp->patch(p).coefsSize();
        }

        // Set shifts
        if (i > 0){
            m_shift_free[i] = m_shift_free[i-1] + m_mappers[i-1].freeSize();
            m_shift_all[i] = m_shift_all[i-1] + m_mappers[i-1].size();
        }
    }

    // FIXIT: move to a seperate method?
    // Save whether boundaries are fixed or tagged
    for (index_t i = 0; i < m_mp->nBoundary(); i++){
        // Get boundary local indices
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);


        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            m_isBoundaryFixed(d,i) = true;
            m_isBoundaryTagged(d,i) = true;
            for (index_t k = 0; k < boundaryDofs.size(); k++)
            {
                if ( m_mappers[d].is_free(boundaryDofs[k],ps.patch) )
                {
                    // If one of the cps on the interface is free,
                    //   we predict that the entire one is free
                    m_isBoundaryFixed(d,i) = false;
                }

                if ( !m_mappers[d].is_tagged(boundaryDofs[k],ps.patch) )
                {
                    // If one of the cps on the interface is not tagged,
                    //  we predict that the entire one is not tagged
                    m_isBoundaryTagged(d,i) = false;
                }

            }
        }
    }

    // Save whether interfaced are fixed or tagged
    for (index_t i = 0; i < m_mp->nInterfaces(); i++){
        // Get boundary local indices
        boundaryInterface interface = m_mp->bInterface(i);
        patchSide ps = interface.first();

        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            m_isInterfaceFixed(d,i) = true;
            m_isInterfaceTagged(d,i) = true;

            for (index_t k = 0; k < boundaryDofs.size(); k++)
            {
                if ( m_mappers[d].is_free(boundaryDofs[k],ps.patch) )
                {
                    // If one of the cps on the interface is free,
                    //   we predict that the entire one is free
                    m_isInterfaceFixed(d,i) = false;
                }

                if ( !m_mappers[d].is_tagged(boundaryDofs[k],ps.patch) )
                {
                    // If one of the cps on the interface is not tagged,
                    //  we predict that the entire one is not tagged
                    m_isInterfaceTagged(d,i) = false;
                }

            }
        }
    }
}

gsVector<> gsParamMethod::getTagged()
{
    gsVector<> out(n_tagged);
    index_t ind = 0;

    for(index_t d = 0; d < m_mp->targetDim(); ++d)
    {
        const std::vector< index_t > t = m_mappers[d].getTagged(); // Get tagged from mapper
        for (std::vector< index_t >::const_iterator it = t.begin(); it != t.end(); ++it)
        {
            // *it returns the global index from the mapper

            // Get local DoFs
            std::vector<std::pair<index_t,index_t> > result;
            m_mappers[d].preImage(*it, result);

            // Get local index and patch
            index_t p = result[0].first;
            index_t i = result[0].second;


            // get the cps here
            out[ind] = m_mp->patch(p).coef(i,d);
            ind++;
        }
    }
    return out;
}

void gsParamMethod::updateTagged(gsVector<> x) const
{
    index_t ind = 0;
    for(index_t d = 0; d < m_mp->targetDim(); ++d)
    {
        const std::vector< index_t > t = m_mappers[d].getTagged(); // Get tagged from mapper
        for (std::vector< index_t >::const_iterator it = t.begin(); it != t.end(); ++it)
        {
            // *it returns the global index from the mapper

            // Get local DoFs
            std::vector<std::pair<index_t,index_t> > result;
            m_mappers[d].preImage(*it, result);

            for(std::vector<std::pair<index_t,index_t> >::const_iterator it = result.begin(); it != result.end(); ++it)
            {
                // Get local index and patch
                index_t p = it->first;
                index_t i = it->second;

                // get the cps here
                m_mp->patch(p).coef(i,d) = x[ind];
            }
            ind++;
        }
    }
}

void gsParamMethod::setupDefaultMappers()
{
    // Get mappers from multibasis with interfaces glued
    gsMultiBasis<> geoBasis(*m_mp);
    for (index_t d = 0; d < m_mp->targetDim(); d++)
    {
        geoBasis.getMapper(iFace::glue,m_mappers[d],false); // False means that we do not finalize
    }

    // Fix coordinates of boundaries
    for (index_t i = 0; i < m_mp->nBoundary(); i++){
        // Get boundary local indices
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            m_mappers[d].markBoundary(ps.patch,boundaryDofs);
        }
    }

    // Finalize mappers
    for (index_t d = 0; d < m_mp->targetDim(); d++)
    {
        m_mappers[d].finalize();
    }

    // Tag coordinates of boundaries
    for (index_t i = 0; i < m_mp->nBoundary(); i++){
        // Get boundary local indices
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            // Tag the controlpoint
            for (index_t j = 0; j < boundaryDofs.size(); j ++){
                m_mappers[d].markTagged(boundaryDofs[j],ps.patch);
            }
        }
    }

    gsInfo << "n controlpoints: " << m_mappers[0].mapSize() << "\n";


    // Call setMappers to calculate n_free, n_tagged etc.
    setMappers(m_mappers);

}

gsVector<> gsParamMethod::getFree() const
{
    // Extract the global design variables
    gsVector<> out(n_free);

    // For all patches
    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        for (index_t i = 0; i < m_mp->patch(p).coefsSize(); i++){
            for (index_t d = 0; d < m_mp->targetDim(); d++){
                index_t ii = m_mappers[d].index(i,p); // Global index
                if (m_mappers[d].is_free_index(ii))
                {
                    out[ii + m_shift_free[d]] = m_mp->patch(p).coef(i,d);
                }
            }
        }
    }

    return out;
}

void gsParamMethod::updateFree(gsVector<> des) const
{
    // For all patches
    for (index_t p = 0; p < m_mp->nBoxes(); p++){
        gsMatrix<> coefs = m_mp->patch(p).coefs(); // Get control points

        for (index_t i = 0; i < coefs.rows(); i++){
            for (index_t d = 0; d < m_mp->targetDim(); d++){
                index_t ii = m_mappers[d].index(i,p); // Global index
                if (m_mappers[d].is_free_index(ii))
                {
                    coefs(i,d) = des[ii + m_shift_free[d]];
                }
            }
        }
        m_mp->patch(p).setCoefs(coefs);
    }
}

void gsParamMethod::updateFreeAndTagged(gsVector<> des, gsVector<> x) const
{
    // Update design variables, i.e. free Dofs
    updateFree(des);

    // Update tagged
    updateTagged(x);

}

gsVector<> gsParamMethod::getControlPoints() const
{

    gsVector<> out(n_cps);
    for(index_t p = 0; p < m_mp->nBoxes(); p++){
        for(index_t i = 0; i < m_mp->patch(p).coefsSize(); i++){
            for(index_t d = 0; d < m_mp->targetDim(); d++){
                index_t ii = m_mappers[d].index(i,p) + m_shift_all[d];
                out[ii] = m_mp->patch(p).coef(i,d);
            }
        }
    }
    return out;
}

gsVector<> gsParamMethod::getControlPoints(gsVector<> des) const
{

    gsVector<> out(n_cps);
    for(index_t p = 0; p < m_mp->nBoxes(); p++){
        for(index_t i = 0; i < m_mp->patch(p).coefsSize(); i++){
            for(index_t d = 0; d < m_mp->targetDim(); d++){
                index_t ii = m_mappers[d].index(i,p);
                index_t gl = ii + m_shift_all[d];

                if (m_mappers[d].is_free_index(ii)){
                    out[gl] = des(ii + m_shift_free[d]);
                } else {
                    out[gl] = m_mp->patch(p).coef(i,d);
                }
            }
        }
    }
    return out;
}

gsVector<> gsParamMethod::getFlat() const
{

    gsVector<> out(n_flat);
    index_t j = 0;
    for ( index_t d = 0; d < m_mp->targetDim(); d++)
    {
        for ( index_t p = 0; p < m_mp->nBoxes(); p++ )
        {
            for (index_t i = 0; i < m_mp->patch(p).coefsSize(); i++)
            {
                out[j] = m_mp->patch(p).coef(i,d);
                j++;
            }
        }
    }
    return out;
}

void gsParamMethod::updateFlat(gsVector<> flat) const
{
    index_t j = 0;
    for ( index_t d = 0; d < m_mp->targetDim(); d++)
    {
        for ( index_t p = 0; p < m_mp->nBoxes(); p++ )
        {
            for (index_t i = 0; i < m_mp->patch(p).coefsSize(); i++)
            {
                m_mp->patch(p).coef(i,d) = flat[j];
                j++;
            }
        }
    }
}

void gsParamMethod::updateControlPoints(gsVector<> cps)
{
    for(index_t p = 0; p < m_mp->nBoxes(); p++){
        for(index_t i = 0; i < m_mp->patch(i).coefsSize(); i++){
            for(index_t d = 0; d < m_mp->targetDim(); d++){
                index_t  ii = m_mappers[d].index(i,p) + m_shift_all[d];
                m_mp->patch(p).coef(i,d) = cps[ii];
            }
        }
    }

}

// FIXIT: looses sparse structure, implement to keep sparsity structure, the
// structure can be predicted by adding two matrices with 1's at nonzeros
gsIpOptSparseMatrix gsParamMethod::mapMatrix(gsDofMapper mapper_in, gsIpOptSparseMatrix M) const
{
    gsMatrix<> mat_in = M.asDense();

    gsMatrix<> mat_out = mapMatrix(mapper_in,mat_in);

    gsIpOptSparseMatrix out(mat_out,-1);

    return out;

};

gsMatrix<> gsParamMethod::mapMatrix(gsDofMapper mapper_in, gsMatrix<> mat_in) const
{
    // Set shifts of input mapper
    gsVector<> mapper_in_shifts;
    mapper_in_shifts.setZero(m_mp->targetDim());

    for(index_t d = 1; d < m_mp->targetDim(); d++){
        mapper_in_shifts[d] = mapper_in_shifts[d-1] + mapper_in.freeSize();
    }

    // FIXIT: take this information as input instead..
    bool row = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.rows();
    bool col = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.cols();
    gsMatrix<> mat_out;

    if(row){
        mat_out.setZero(n_free,mat_in.cols());
    } else if (col) {
        mat_out.setZero(mat_in.rows(),n_free);
    } else {
        GISMO_ERROR("Wrong input size in mapMatrix..\n");
    }

    for(index_t d = 0; d < m_mp->targetDim(); d++){
        // Iterate through free indices
        for (index_t ii = 0; ii < m_mappers[d].freeSize(); ii++){
            // Get a local index
            std::vector<std::pair<index_t,index_t> > result;
            m_mappers[d].preImage(ii, result);

            for(std::vector<std::pair<index_t,index_t>>::iterator it=result.begin(); it != result.end(); ++it)
            {
                // Get local index and patch
                index_t p = it->first;
                index_t i = it->second;

                // Convert to global to find the right column
                index_t ii2 = mapper_in.index(i,p) + mapper_in_shifts[d];

                if (row){ // If right no colms
                    mat_out.row(ii + m_shift_free[d]) += mat_in.row(ii2);
                } else if (col){
                    mat_out.col(ii + m_shift_free[d]) += mat_in.col(ii2);
                }
            }
        }
    }
    return mat_out;

};

gsMatrix<> gsParamMethod::mapMatrixToTagged(gsDofMapper mapper_in, gsMatrix<> mat_in) const
{
    // Set shifts of input mapper
    gsVector<> mapper_in_shifts, tagged_shift;
    mapper_in_shifts.setZero(m_mp->targetDim());
    tagged_shift.setZero(m_mp->targetDim());

    for(index_t d = 1; d < m_mp->targetDim(); d++){
        mapper_in_shifts[d] = mapper_in_shifts[d-1] + mapper_in.freeSize();
        tagged_shift[d] = tagged_shift[d-1] + m_mappers[d].taggedSize();
    }

    // FIXIT: take this information as input instead..
    bool row = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.rows();
    bool col = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.cols();
    gsMatrix<> mat_out;

    if(row){
        mat_out.setZero(n_tagged,mat_in.cols());
    } else if (col) {
        mat_out.setZero(mat_in.rows(),n_tagged);
    } else {
        GISMO_ERROR("Wrong input size in mapMatrix..\n");
    }

    for(index_t d = 0; d < m_mp->targetDim(); d++){
        // Iterate through free indices
        for (index_t t = 0; t < m_mappers[d].taggedSize(); t++){
            // Get global index
            index_t ii = m_mappers[d].getTagged()[t];
            // Get a local index
            std::vector<std::pair<index_t,index_t> > result;
            m_mappers[d].preImage(ii, result);

            for(std::vector<std::pair<index_t,index_t>>::iterator it=result.begin(); it != result.end(); ++it)
            {
                // Get local index and patch
                index_t p = it->first;
                index_t i = it->second;

                // Convert to global to find the right column
                index_t ii2 = mapper_in.index(i,p) + mapper_in_shifts[d];

                if (row){ // If right no colms
                    mat_out.row(t + tagged_shift[d]) += mat_in.row(ii2);
                } else if (col){
                    mat_out.col(t + tagged_shift[d]) += mat_in.col(ii2);
                }
            }
        }
    }
    return mat_out;

};

// FIXIT: Implement refExtension functionally.. ask Angelos
void gsParamMethod::refineElements(const std::vector<bool> & elMarked)
{
    gsMultiBasis<> basis(*m_mp);
    const int dim = basis.dim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    int numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<> refBoxes;

    for (unsigned pn=0; pn < basis.nBases(); ++pn )// for all patches
    {
        // Get number of elements to be refined on this patch
        const int numEl = basis[pn].numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  std::bind2nd(std::equal_to<bool>(), true) );
        poffset += numEl;
        refBoxes.resize(dim, 2*numMarked);
        //gsDebugVar(numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<>::domainIter domIt = basis.basis(pn).makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount++ ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                refBoxes.col(2*numMarked  ) =
                        refBoxes.col(2*numMarked+1) = domIt->centerPoint();

                // Advance marked cells counter
                numMarked++;
            }
        }
        // Refine all of the found refBoxes in this patch
        // FIXIT: make dimension independet.. Find way to replace 2 with dim.
        std::vector<gsSortedVector<unsigned> > OX = dynamic_cast<gsHTensorBasis<2,real_t>&> (m_mp->patch(pn).basis()).getXmatrix();
        m_mp->patch(pn).basis().refine( refBoxes );
        gsSparseMatrix<> transf;
        dynamic_cast<gsHTensorBasis<2,real_t>&> (m_mp->patch(pn).basis()).transfer(OX, transf);
        //gsDebug<<"tranf orig:\n"<<transf<<std::endl;
        m_mp->patch(pn).coefs() = transf*m_mp->patch(pn).coefs();

    }

    m_mp->repairInterfaces();
}

void gsParamMethod::recreateMappers()
{
    // Get mappers from multibasis with interfaces glued
    gsMultiBasis<> geoBasis(*m_mp);
    for (index_t d = 0; d < m_mp->targetDim(); d++)
    {
        geoBasis.getMapper(iFace::glue,m_mappers[d],false); // False means that we do not finalize
    }

    // Check if we should fix or tag some boundaries
    for (index_t i = 0; i < m_mp->nBoundary(); i++){
        // Get boundary local indices
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            if (m_isBoundaryFixed(d,i))
                m_mappers[d].markBoundary(ps.patch,boundaryDofs);
        }
    }

    // Check if we should fix or tag some interfaces
    for (index_t i = 0; i < m_mp->nInterfaces(); i++){
        // Get boundary local indices
        boundaryInterface interface = m_mp->bInterface(i);
        patchSide ps = interface.first();

        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            if (m_isInterfaceFixed(d,i))
                m_mappers[d].markBoundary(ps.patch,boundaryDofs);
        }
    }

    // Finalize mappers
    for (index_t d = 0; d < m_mp->targetDim(); d++)
    {
        m_mappers[d].finalize();
    }

    // Tag coordinates of boundaries
    for (index_t i = 0; i < m_mp->nBoundary(); i++){
        // Get boundary local indices
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            // Tag the controlpoint
            if (m_isBoundaryTagged(d,i))
            {
                for (index_t j = 0; j < boundaryDofs.size(); j ++){
                    m_mappers[d].markTagged(boundaryDofs[j],ps.patch);
                }
            }
        }
    }

    // Tag coordinates of interfaces
    for (index_t i = 0; i < m_mp->nInterfaces(); i++){
        // Get boundary local indices
        boundaryInterface interface = m_mp->bInterface(i);
        patchSide ps = interface.first();

        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            // Tag the controlpoint
            if (m_isInterfaceTagged(d,i))
            {
                for (index_t j = 0; j < boundaryDofs.size(); j ++){
                    m_mappers[d].markTagged(boundaryDofs[j],ps.patch);
                }
            }
        }
    }

    gsInfo << "n controlpoints: " << m_mappers[0].mapSize() << "\n";


    // Call setMappers to calculate n_free, n_tagged etc.
    setMappers(m_mappers);

}
