#include <gismo.h>
#include <fstream>
#include "gsParamMethod.h"
using namespace gismo;

gsParamMethod::gsParamMethod(gsMultiPatch<>* mpin,std::vector< gsDofMapper > mappers):
    m_mp(mpin)
{
    m_mappers = mappers;

    m_shift_free.setZero(m_mp->targetDim());
    m_shift_all.setZero(m_mp->targetDim());

    n_free = 0;
    n_cps = 0;
    n_tagged = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_free += m_mappers[i].freeSize();
        n_cps += m_mappers[i].size();
        n_tagged += m_mappers[i].taggedSize();

        // Set shifts
        if (i > 0){
            m_shift_free[i] = m_shift_free[i-1] + m_mappers[i-1].freeSize();
            m_shift_all[i] = m_shift_all[i-1] + m_mappers[i-1].size();
        }
    }
}

gsParamMethod::gsParamMethod(gsMultiPatch<>* mpin): m_mp(mpin), m_mappers(m_mp->targetDim())
{
    setupMapper(); // FIXIT should be changed to general setup rather than hardcoded for antenna
    m_shift_free.setZero(m_mp->targetDim());
    m_shift_all.setZero(m_mp->targetDim());

    n_free = 0;
    n_cps = 0;
    n_tagged = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_free += m_mappers[i].freeSize();
        n_cps += m_mappers[i].size();
        n_tagged += m_mappers[i].taggedSize();

        // Set shifts
        if (i > 0){
            m_shift_free[i] = m_shift_free[i-1] + m_mappers[i-1].freeSize();
            m_shift_all[i] = m_shift_all[i-1] + m_mappers[i-1].size();
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

void gsParamMethod::setupMapper()
{
  //FIXIT: Patch 3 is hardcoded in paraOptProblem

  // Get mappers from multibasis with interfaces glued
  gsMultiBasis<> geoBasis(*m_mp);
  geoBasis.getMapper(iFace::glue,m_mappers[0],false); // False means that we do not finalize
  geoBasis.getMapper(iFace::glue,m_mappers[1],false); // False means that we do not finalize

  // Fix y coordinate for all boundaries
  for (index_t i = 0; i < m_mp->nBoundary(); i ++){
    patchSide ps = m_mp->boundaries()[i];
    gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
    // gsInfo << ps << " is fixed in y direction \n\n";
    // gsInfo << "boundaryDofs on patch " << ps.patch << "\n"<< boundaryDofs << "\n\n";
    m_mappers[1].markBoundary(ps.patch,boundaryDofs);

    // FIXIT it is hardcoded that only xcoordinates of top bnd is fixed
    if (m_mp->nBoxes() > 3){
        if (ps.patch == 4){
            m_mappers[0].markBoundary(ps.patch,boundaryDofs);
        }
    } else {
        // If there is less than 3 patches then all boundaries is fixed
        m_mappers[0].markBoundary(ps.patch,boundaryDofs);
    }
  }

  // Fix boundary of the fixedPatch
  gsMultiPatch<> fixedGeom(m_mp->patch(fixedPatch));

  for(index_t i = 0; i < fixedGeom.nBoundary(); i++){
    patchSide ps = fixedGeom.boundaries()[i];
    gsVector<unsigned> boundaryDofs = m_mp->basis(fixedPatch).boundary(ps);
    m_mappers[0].markBoundary(fixedPatch,boundaryDofs); // Mark xcoord of boundaries of fixed patch
    m_mappers[1].markBoundary(fixedPatch,boundaryDofs); // Mark ycoord of boundaries of
  }

  // Finalize mappers
  m_mappers[0].finalize();
  m_mappers[1].finalize();

  n_controlpoints = m_mappers[0].mapSize();
  gsInfo << "n cps: " << n_controlpoints << "\n";

  // Tag fixed bnds
  for(index_t i = 0; i < fixedGeom.nBoundary(); i++){
    patchSide ps = fixedGeom.boundaries()[i];
    gsVector<unsigned> boundaryDofs = m_mp->basis(fixedPatch).boundary(ps);
    for (index_t j = 0; j < boundaryDofs.size(); j ++){
        m_mappers[0].markTagged(boundaryDofs[j],fixedPatch); // Mark xcoord of boundaries of fixed patch
        m_mappers[1].markTagged(boundaryDofs[j],fixedPatch); // Mark ycoord of boundaries of

    }
  }

  // count number of free dofs
  n_free = 0;
  for(index_t i = 0; i < m_mp->targetDim(); i++){
      n_free += m_mappers[i].freeSize();
  }
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
        for(index_t i = 0; i < m_mp->patch(i).coefsSize(); i++){
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
