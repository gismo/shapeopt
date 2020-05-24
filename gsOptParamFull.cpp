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
        index_t k2 = boundaryDofs1[i];

        for( index_t d = 0; d < m_mp->targetDim(); d++)
        {
            result += pow( m_mp->patch(ps1.patch).coef(k1,d) - m_mp->patch(ps2.patch).coef(k2,d) , 2);
        }
    }

    return result;
    // Why does it return 0?


}

real_t gsOptParamFull::evalObj() const {
    gsVector<> diff = (m_tagged_goal - m_paramMethod->getTagged());

    real_t result = 0.5*diff.squaredNorm();
    gsDebugVar(result);

    for(index_t i = 0; i < m_mp->nInterfaces(); i++)
    {
        boundaryInterface interface = m_mp->bInterface(i);

        result += 0.5*getObjFromInterface(interface, m_mp);
        gsDebugVar(result);
    }

    gsDebugVar(result);
    gsInfo << "DONE\n" << std::flush;
    return result;
}

gsVector<> gsOptParamFull::gradAll() const{
    gsVector<> out;
    out.setZero(n_flat);

    index_t iter = 0;
    for (index_t d = 0; d < m_mp->targetDim(); d++){
        for (index_t k = 0; k < m_mp->nBoxes(); k++){
            for (index_t i = 0; i < m_mp->patch(k).coefsSize(); i++){
                if (m_mappers[d].is_tagged(i,k))
                {
                    index_t ii = m_mappers[d].index(i,k);
                    std::vector< std::pair < index_t, index_t > > res;

                    m_mappers[d].preImage(ii,res);

                    real_t size = res.size();
                    // real_t size = 1;

                    out[iter] = 1.0/size * (m_mp->patch(k).coef(i,d) - m_mp_goal->patch(k).coef(i,d));
                }
                iter++;
            }
        }
    }


    return out;

};

void gsOptParamFull::setupMappers()
{
    // Get mappers from multibasis with interfaces glued
    gsMultiBasis<> geoBasis(*m_mp);
    for (index_t d = 0; d < m_mp->targetDim(); d++)
    {
        gsDofMapper tmp(geoBasis);
        m_mappers[d] = tmp;
        m_mappers[d].finalize();
    }


    // count number of free dofs
    n_free = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_free += m_mappers[i].freeSize();
    }

}

