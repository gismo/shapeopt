#include <gismo.h>
#include <math.h>
#include "gsOptParamFull.h"

real_t getObjFromInterface(boundaryInterface interface, gsMultiPatch<>::Ptr m_mp)
{
    real_t result = 0;
    
    patchSide ps1 = interface.first();
    patchSide ps2 = interface.second();

    gsVector<unsigned> boundaryDofs1 = m_mp->basis(ps1.patch).boundary(ps1);
    gsVector<unsigned> boundaryDofs2 = m_mp->basis(ps2.patch).boundary(ps2);

    GISMO_ASSERT(boundaryDofs1.size() == boundaryDofs2.size(), "Interfaces dont match!");

    for( index_t i = 0; i < boundaryDofs1.size(); i++)
    {
        index_t k1 = boundaryDofs1[i];
        index_t k2 = boundaryDofs2[i];

        for( index_t d = 0; d < m_mp->targetDim(); d++)
        {
            result += 0.5*pow( m_mp->patch(ps1.patch).coef(k1,d) - m_mp->patch(ps2.patch).coef(k2,d) , 2);
        }
    }

    return result;
    // Why does it return 0?


}

gsVector<> gsOptParamFull::getGradAllFromInterface(boundaryInterface interface) const
{
    gsVector<> result;
    result.setZero(n_flat);
    
    patchSide ps1 = interface.first();
    patchSide ps2 = interface.second();

    gsVector<unsigned> boundaryDofs1 = m_mp->basis(ps1.patch).boundary(ps1);
    gsVector<unsigned> boundaryDofs2 = m_mp->basis(ps2.patch).boundary(ps2);

    GISMO_ASSERT(boundaryDofs1.size() == boundaryDofs2.size(), "Interfaces dont match!");

    for( index_t i = 0; i < boundaryDofs1.size(); i++)
    {
        index_t k1 = boundaryDofs1[i];
        index_t k2 = boundaryDofs2[i];

        for( index_t d = 0; d < m_mp->targetDim(); d++)
        {

            index_t ii = k1 + patchShift[ps1.patch] + m_paramMethod->m_shift_flat[d];
            result[ii] += m_mp->patch(ps1.patch).coef(k1,d) - m_mp->patch(ps2.patch).coef(k2,d);

            index_t jj = k2 + patchShift[ps2.patch] + m_paramMethod->m_shift_flat[d];
            result[jj] += -(m_mp->patch(ps1.patch).coef(k1,d) - m_mp->patch(ps2.patch).coef(k2,d));
        }
    }

    return result;
    // Why does it return 0?


}

real_t gsOptParamFull::evalObj() const {
    real_t result = 0;

    for (index_t d = 0; d < m_mp->targetDim(); d++){
        for (index_t k = 0; k < m_mp->nBoxes(); k++){
            for (index_t i = 0; i < m_mp->patch(k).coefsSize(); i++){
                if (m_mappers[d].is_tagged(i,k))
                {
                    result += 0.5*pow(m_mp->patch(k).coef(i,d) - m_mp_goal->patch(k).coef(i,d), 2);
                }
            }
        }
    }
    //gsDebugVar(result);

    for(index_t i = 0; i < m_mp->nInterfaces(); i++)
    {
        boundaryInterface interface = m_mp->bInterface(i);

        result += getObjFromInterface(interface, m_mp);
    }

    //gsDebugVar(result);
    //gsInfo << "DONE\n" << std::flush;
    return result;
}

gsVector<> gsOptParamFull::gradAll() const{
    gsInfo << "Grad All Opt Param Full\n";
    gsVector<> out;
    out.setZero(n_flat);

    index_t iter = 0;
    for (index_t d = 0; d < m_mp->targetDim(); d++){
        for (index_t k = 0; k < m_mp->nBoxes(); k++){
            for (index_t i = 0; i < m_mp->patch(k).coefsSize(); i++){
                if (m_mappers[d].is_tagged(i,k))
                {
                    out[iter] =  (m_mp->patch(k).coef(i,d) - m_mp_goal->patch(k).coef(i,d));
                }
                iter++;
            }
        }
    }

    for(index_t i = 0; i < m_mp->nInterfaces(); i++)
    {
        boundaryInterface interface = m_mp->bInterface(i);

        out += getGradAllFromInterface(interface);
    }


    return out;

};

void gsOptParamFull::setupMappers()
{
    m_mappers_old = m_mappers;
    // Get mappers from multibasis with interfaces glued
    gsMultiBasis<> geoBasis(*m_mp);
    for (index_t d = 0; d < m_mp->targetDim(); d++)
    {
        gsDofMapper tmp(geoBasis);
        m_mappers[d] = tmp;
        m_mappers[d].finalize();
    }

    // Tag again
    for (index_t d = 0; d < m_mp->targetDim(); d++)
    {
        for (index_t k = 0; k < m_mp->nBoxes(); k++)
        {
            for (index_t i = 0; i < m_mp->patch(k).coefsSize(); i++)
            {
               if (m_mappers_old[d].is_tagged(i,k))
               {
                   m_mappers[d].markTagged(i,k);
               }
            }
        }
    }


    // count number of free dofs and patchShift
    n_free = 0;
    n_tagged = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_free += m_mappers[i].freeSize();
        n_tagged += m_mappers[i].taggedSize();
    }

    patchShift.setZero(m_mp->nBoxes());
    index_t tmpflat = 0;
    for(index_t k = 0; k < m_mp->nBoxes(); k++){
        patchShift[k] = tmpflat;
        tmpflat += m_mp->patch(k).coefsSize();
    }

    m_desLowerBounds.setConstant(n_tagged, -100);
    m_desUpperBounds.setConstant(n_tagged, 100);

    m_paramMethod->setMappers(m_mappers);

}
