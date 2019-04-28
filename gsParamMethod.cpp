#include <gismo.h>
#include <fstream>
#include "gsParamMethod.h"
using namespace gismo;

gsParamMethod::gsParamMethod(std::vector< gsDofMapper > mappers)
{
    m_mappers = mappers;
    if (m_mappers.size() > 1){
        GISMO_ASSERT(m_mappers[1].firstIndex() == m_mappers[0].freeSize(),
            "the second mapper must be shifted with the freesize of the first mapper");
    }
    if (m_mappers.size() > 2){
        GISMO_ASSERT(m_mappers[2].firstIndex() == m_mappers[1].freeSize() + m_mappers[0].freeSize(),
            "the third mapper must be shifted with the freesize of the first two mappers");
    }
    n_freeCps = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_freeCps += m_mappers[i].freeSize();
    }
}

gsParamMethod::gsParamMethod(gsMultiPatch<>* mpin): m_mp(mpin), m_mappers(m_mp->targetDim())
{
  setupMapper(); // FIXIT should be changed to general setup rather than hardcoded for antenna
  n_freeCps = 0;
  for(index_t i = 0; i < m_mp->targetDim(); i++){
      n_freeCps += m_mappers[i].freeSize();
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

  // Finalize first mapper
  m_mappers[0].finalize();

  // Shift the y mapper and finalize:
  m_mappers[1].setShift(m_mappers[0].freeSize());
  m_mappers[1].finalize();

  n_controlpoints = m_mappers[0].mapSize();
  gsInfo << "n cps: " << n_controlpoints << "\n";

  n_freeCps = 0;
  for(index_t i = 0; i < m_mp->targetDim(); i++){
      n_freeCps += m_mappers[i].freeSize();
  }
}

gsVector<> gsParamMethod::getDesignVariables() const
{
    // Extract the global design variables
    gsVector<> out(n_freeCps);

    // For all patches
    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    { // For each patch
        gsMatrix<> coefs = m_mp->patch(p).coefs(); // Get control points

        for (index_t i = 0; i < coefs.rows(); i++){
            for (index_t d = 0; d < m_mp->targetDim(); d++){
                index_t ii = m_mappers[d].index(i,p); // Global index
                if (m_mappers[d].is_free_index(ii))
                out[ii] = coefs(i,d);
            }
        }
    }

    return out;
}

void gsParamMethod::updateDesignVariables(gsVector<> des) const
{
    // For all patches
    for (index_t p = 0; p < m_mp->nBoxes(); p++){
        gsMatrix<> coefs = m_mp->patch(p).coefs(); // Get control points

        for (index_t i = 0; i < coefs.rows(); i++){
            for (index_t d = 0; d < m_mp->targetDim(); d++){
                index_t ii = m_mappers[d].index(i,p); // Global index
                if (m_mappers[d].is_free_index(ii))
                coefs(i,d) = des[ii];
            }
        }
        m_mp->patch(p).setCoefs(coefs);
    }
}

// FIXIT: make dimension independent?
gsVector<> gsParamMethod::getControlPoints() const {

    gsVector<> cx(n_controlpoints);
    gsVector<> cy(n_controlpoints);
    index_t j = 0;
    for(index_t i = 0; i < m_mp->nBoxes(); i++){
        for(index_t k = 0; k < m_mp->patch(i).coefsSize(); k++){
            cx[j] = m_mp->patch(i).coef(k,0);
            cy[j] = m_mp->patch(i).coef(k,1);
            j++;
        }
    }

    gsVector<> out(2*n_controlpoints);
    out << cx,cy;

    return out;
}

// FIXIT: make dimension independent?
gsVector<> gsParamMethod::getControlPoints(gsVector<> des) const {
    gsVector<> cx(n_controlpoints);
    gsVector<> cy(n_controlpoints);
    index_t j = 0;
    for(index_t p = 0; p < m_mp->nBoxes(); p++){
        for(index_t i = 0; i < m_mp->patch(p).coefsSize(); i++){
            index_t iix = m_mappers[0].index(i,p);

            if (m_mappers[0].is_free_index(iix)){
                cx[j] = des[iix];
            } else {
                cx[j] = m_mp->patch(p).coef(i,0);
            }

            index_t iiy = m_mappers[1].index(i,p);

            if (m_mappers[1].is_free_index(iiy)){
                cy[j] = des[iiy];
            } else {
                cy[j] = m_mp->patch(p).coef(i,1);
            }
            j++;
        }
    }

    gsVector<> out(2*n_controlpoints);
    out << cx,cy;

    return out;
}
