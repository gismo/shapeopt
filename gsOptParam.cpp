#include <gismo.h>
#include "gsOptParam.h"

// Implement
gsOptParam::gsOptParam(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsMultiPatch<>> mp_goal, memory::shared_ptr<gsShapeOptLog> slog, index_t param, bool use_Lagrangian):
    gsShapeOptProblem(mp,slog),
    m_mp_goal(mp_goal),
    m_pM_goal(mp_goal)
{
    m_tagged_goal = m_pM_goal.getTagged();

    std::string name = "goal";
    m_log->saveVec(m_pM_goal.getFlat(),name);

    setupMappers();
    // Allocate parametrization method
    if (param == 0){
        m_paramMethod = memory::make_shared(new gsSpringMethod(m_mp,m_mappers));
        m_paramMethod->computeMap();
    } else if (param == 1) {
        gsModLiao::Ptr opt_param = memory::make_shared(new gsModLiao(m_mp,m_mappers,true));
        opt_param->setQuad(m_quA,m_quB);
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 2) {
        gsWinslow::Ptr opt_param = memory::make_shared(new gsWinslow(m_mp,m_mappers,true));
        opt_param->setQuad(m_quA,m_quB);// FIXIT: this and the next two statements can be moved out of if state ment if param =! 0 if ()...
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 3) {
        gsLiao::Ptr opt_param = memory::make_shared(new gsLiao(m_mp,m_mappers,true));
        opt_param->setQuad(m_quA,m_quB);
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 4) {
        gsHarmonic::Ptr opt_param = memory::make_shared(new gsHarmonic(m_mp,m_mappers,true));
        opt_param->setQuad(m_quA,m_quB);
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 5) {
        m_paramMethod = memory::make_shared(new gsWinslowWithDeriv(m_mp,m_mappers,false,false,true,0));
        // m_paramMethod->setQuad(m_quA,m_quB); FIXIT: implement, maybe with dynamic cast
    } else if (param == 6) {
        gsWinslow::Ptr opt_param = memory::make_shared(new gsWinslow(m_mp,m_mappers,false,false,true,0));
        // opt_param->setQuad(m_quA,m_quB);// FIXIT: this and the next two statements can be moved out of if state ment if param =! 0 if ()...
        opt_param->setQuad(m_quA,m_quB);
        *m_log << "quA, quB = " << opt_param->m_quA << ", " << opt_param->m_quB << "\n";
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else {
        GISMO_ERROR("Param value is not 0 to 6... wrong input..\n");
    }

    // Copy some data from m_paramMethod for ease of use
    n_free = m_paramMethod->n_free;
    n_flat = m_paramMethod->n_flat;
    n_tagged = m_paramMethod->n_tagged;
    n_cps = m_paramMethod->n_cps;

    name = "init";
    m_log->saveVec(m_paramMethod->getTagged(),name);

    // Also calls setupDesignBounds
    setupOptParameters();

}

real_t gsOptParam::evalObj() const {
    gsVector<> diff = (m_tagged_goal - m_paramMethod->getTagged());
    return 0.5*diff.squaredNorm();
}

gsVector<> gsOptParam::gradObj() const{
    return m_paramMethod->getTagged() - m_tagged_goal;
};

gsVector<> gsOptParam::gradAll() const{
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

void gsOptParam::setupDesignBounds()
{
    m_desLowerBounds.setZero(n_tagged);
    // X coordinates
    gsVector<> xLowerBounds;
    xLowerBounds.setConstant(m_mappers[0].taggedSize(),-100);
    m_desLowerBounds.segment(0,m_mappers[0].taggedSize()) = xLowerBounds;

    // Y coordinates
    gsVector<> yLowerBounds;
    yLowerBounds.setConstant(m_mappers[1].taggedSize(),-100);
    m_desLowerBounds.segment(m_mappers[0].taggedSize(),m_mappers[1].taggedSize()) = yLowerBounds;

    if (m_mp->targetDim() == 3)
    {
        // Z coordinates
        gsVector<> zLowerBounds;
        zLowerBounds.setConstant(m_mappers[2].taggedSize(),-100);
        m_desLowerBounds.segment(m_mappers[1].taggedSize() + m_mappers[0].taggedSize(),m_mappers[2].taggedSize()) = zLowerBounds;
    }

    m_desUpperBounds.setZero(n_tagged);
    // X coordinates
    gsVector<> xUpperBounds;
    xUpperBounds.setConstant(m_mappers[0].taggedSize(),100);
    m_desUpperBounds.segment(0,m_mappers[0].taggedSize()) = -xLowerBounds;

    // Y coordinates
    gsVector<> yUpperBounds;
    yUpperBounds.setConstant(m_mappers[1].taggedSize(),100);
    m_desUpperBounds.segment(m_mappers[0].taggedSize(),m_mappers[1].taggedSize()) = yUpperBounds;

    if (m_mp->targetDim() == 3)
    {
        // Z coordinates
        gsVector<> zUpperBounds;
        zUpperBounds.setConstant(m_mappers[2].taggedSize(),100);
        m_desUpperBounds.segment(m_mappers[1].taggedSize() + m_mappers[0].taggedSize(),m_mappers[2].taggedSize()) = zUpperBounds;
    }
}

void gsOptParam::setupMappers(){
    m_mappers = m_pM_goal.mappers();
}

gsDofMapper gsOptParam::mapper_grad() const
{
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    gsExprAssembler<>::space u = A.getSpace(dbasis,m_mp->geoDim()); // The gradient is a vector with targetDim values

    A.initSystem();

    return u.mapper();
}
