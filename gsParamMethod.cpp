#include <gismo.h>
#include <fstream>
#include "paraWithMapper.h"
using namespace gismo;

paraWithMapper::paraWithMapper(gsMultiPatch<>* mpin): mp(mpin), m_mappers(mp->targetDim()){
};

void paraWithMapper::setupMapper(){
  //FIXIT: Patch 3 is hardcoded in paraOptProblem

  // Get mappers from multibasis with interfaces glued
  gsMultiBasis<> geoBasis(*mp);
  geoBasis.getMapper(iFace::glue,m_mappers[0],false); // False means that we do not finalize
  geoBasis.getMapper(iFace::glue,m_mappers[1],false); // False means that we do not finalize

  // Fix y coordinate for all boundaries
  for (index_t i = 0; i < mp->nBoundary(); i ++){
    patchSide ps = mp->boundaries()[i];
    gsVector<unsigned> boundaryDofs = mp->basis(ps.patch).boundary(ps);
    // gsInfo << ps << " is fixed in y direction \n\n";
    // gsInfo << "boundaryDofs on patch " << ps.patch << "\n"<< boundaryDofs << "\n\n";
    m_mappers[1].markBoundary(ps.patch,boundaryDofs);

    // FIXIT it is hardcoded that only xcoordinates of top bnd is fixed
    if (mp->nBoxes() > 3){
        if (ps.patch == 4){
            m_mappers[0].markBoundary(ps.patch,boundaryDofs);
        }
    } else {
        // If there is less than 3 patches then all boundaries is fixed
        m_mappers[0].markBoundary(ps.patch,boundaryDofs);
    }
  }

  // Fix boundary of the fixedPatch
  gsMultiPatch<> fixedGeom(mp->patch(fixedPatch));

  for(index_t i = 0; i < fixedGeom.nBoundary(); i++){
    patchSide ps = fixedGeom.boundaries()[i];
    gsVector<unsigned> boundaryDofs = mp->basis(fixedPatch).boundary(ps);
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
  for(index_t i = 0; i < mp->targetDim(); i++){
      n_freeCps += m_mappers[i].freeSize();
  }
}

gsVector<> paraWithMapper::getDesignVariables() const{
  // Extract the global design variables
  gsVector<> out(m_mappers[0].freeSize() + m_mappers[1].freeSize());

  // For all patches
  for (index_t p = 0; p < mp->nBoxes(); p++){
     gsMatrix<> coefs = mp->patch(p).coefs(); // Get control points

     for (index_t i = 0; i < coefs.rows(); i++){
         index_t iix = m_mappers[0].index(i,p); // Global index
         if (m_mappers[0].is_free_index(iix))
            out[iix] = coefs(i,0); // x coordinate

         index_t iiy = m_mappers[1].index(i,p); // Global index
         if (m_mappers[1].is_free_index(iiy))
            out[iiy] = coefs(i,1); // y coordinate
     }
  }

  return out;
}

void paraWithMapper::updateDesignVariables(gsVector<> des) const {
  // For all patches
  for (index_t p = 0; p < mp->nBoxes(); p++){
     gsMatrix<> coefs = mp->patch(p).coefs(); // Get control points

     for (index_t i = 0; i < coefs.rows(); i++){
         index_t iix = m_mappers[0].index(i,p); // Global index
         if (m_mappers[0].is_free_index(iix))
            coefs(i,0) = des[iix]; // x coordinate

         index_t iiy = m_mappers[1].index(i,p); // Global index
         if (m_mappers[1].is_free_index(iiy))
            coefs(i,1) = des[iiy]; // y coordinate
     }
     mp->patch(p).setCoefs(coefs);
  }
}

gsVector<> paraWithMapper::getControlPoints() const {

    gsVector<> cx(n_controlpoints);
    gsVector<> cy(n_controlpoints);
    index_t j = 0;
    for(index_t i = 0; i < mp->nBoxes(); i++){
        for(index_t k = 0; k < mp->patch(i).coefsSize(); k++){
            cx[j] = mp->patch(i).coef(k,0);
            cy[j] = mp->patch(i).coef(k,1);
            j++;
        }
    }

    gsVector<> out(2*n_controlpoints);
    out << cx,cy;

    return out;
}

gsVector<> paraWithMapper::getControlPoints(gsVector<> des) const {
    gsVector<> cx(n_controlpoints);
    gsVector<> cy(n_controlpoints);
    index_t j = 0;
    for(index_t p = 0; p < mp->nBoxes(); p++){
        for(index_t i = 0; i < mp->patch(p).coefsSize(); i++){
            index_t iix = m_mappers[0].index(i,p);

            if (m_mappers[0].is_free_index(iix)){
                cx[j] = des[iix];
            } else {
                cx[j] = mp->patch(p).coef(i,0);
            }

            index_t iiy = m_mappers[1].index(i,p);

            if (m_mappers[1].is_free_index(iiy)){
                cy[j] = des[iiy];
            } else {
                cy[j] = mp->patch(p).coef(i,1);
            }
            j++;
        }
    }

    gsVector<> out(2*n_controlpoints);
    out << cx,cy;

    return out;
}
