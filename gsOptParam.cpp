#include <gismo.h>
#include "gsOptParam.h"

// Implement
gsOptParam::gsOptParam(gsMultiPatch<>* mp, gsMultiPatch<>* mp_goal, gsShapeOptLog* slog, index_t param, bool use_Lagrangian):
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
        m_paramMethod = new gsSpringMethod(m_mp,m_mappers);
        m_paramMethod->computeMap();
    } else if (param == 1) {
        gsModLiao *opt_param = new gsModLiao(m_mp,m_mappers,true);
        // opt_param->setQuad(quA,quB);
        m_paramMethod = new gsAffineOptParamMethod(opt_param, use_Lagrangian);
        m_paramMethod->computeMap();
    } else if (param == 2) {
        gsWinslow *opt_param = new gsWinslow(m_mp,m_mappers,true);
        // opt_param->setQuad(quA,quB);// FIXIT: this and the next two statements can be moved out of if state ment if param =! 0 if ()...
        m_paramMethod = new gsAffineOptParamMethod(opt_param, use_Lagrangian);
        m_paramMethod->computeMap();
    } else if (param == 3) {
        gsLiao *opt_param = new gsLiao(m_mp,m_mappers,true);
        // opt_param->setQuad(quA,quB);
        m_paramMethod = new gsAffineOptParamMethod(opt_param, use_Lagrangian);
        m_paramMethod->computeMap();
    } else if (param == 4) {
        gsHarmonic *opt_param = new gsHarmonic(m_mp,m_mappers,true);
        // opt_param->setQuad(quA,quB);
        m_paramMethod = new gsAffineOptParamMethod(opt_param, use_Lagrangian);
        m_paramMethod->computeMap();
    } else if (param ==5) {
        m_paramMethod = new gsWinslowWithDeriv(m_mp,m_mappers,false,false,true,0);
        // m_paramMethod->setQuad(quA,quB); FIXIT: implement, maybe with dynamic cast
    } else {
        GISMO_ERROR("Param value is not 0 to 5... wrong input..\n");
    }

    // Copy some data from m_paramMethod for ease of use
    n_free = m_paramMethod->n_free;
    n_flat = m_paramMethod->n_flat;
    n_tagged = m_paramMethod->n_tagged;
    n_cps = m_paramMethod->n_cps;

    setupOptParameters();

    m_desLowerBounds.setZero(n_tagged);
    // X coordinates
    gsVector<> xLowerBounds;
    xLowerBounds.setConstant(m_mappers[0].taggedSize(),-100);
    m_desLowerBounds.segment(0,m_mappers[0].taggedSize()) = xLowerBounds;

    // Y coordinates
    gsVector<> yLowerBounds;
    yLowerBounds.setConstant(m_mappers[1].taggedSize(),-100);
    m_desLowerBounds.segment(m_mappers[0].taggedSize(),m_mappers[1].taggedSize()) = yLowerBounds;

    m_desUpperBounds.setZero(n_tagged);
    // X coordinates
    gsVector<> xUpperBounds;
    xUpperBounds.setConstant(m_mappers[0].taggedSize(),100);
    m_desUpperBounds.segment(0,m_mappers[0].taggedSize()) = -xLowerBounds;

    // Y coordinates
    gsVector<> yUpperBounds;
    yUpperBounds.setConstant(m_mappers[1].taggedSize(),100);
    m_desUpperBounds.segment(m_mappers[0].taggedSize(),m_mappers[1].taggedSize()) = yUpperBounds;
}

real_t gsOptParam::evalObj() const {
    gsVector<> diff = (m_tagged_goal - m_paramMethod->getTagged());
    return 0.5*diff.squaredNorm();
}

gsVector<> gsOptParam::gradObj() const{
    return m_paramMethod->getTagged() - m_tagged_goal;
};

void gsOptParam::setupMappers(){
    m_mappers = m_pM_goal.mappers();
}
