#include <gismo.h>
#include "gsAffineOptParamMethod.h"
using namespace gismo;

gsAffineOptParamMethod::gsAffineOptParamMethod(gsOptParamMethod* optParamMethod):
    gsAffineParamMethod(optParamMethod->mp(),optParamMethod->mappers()),m_optParamMethod(optParamMethod)
{
    reset(); //Setup the problem with the parametrization hold on m_pM->m_mp as reference
};

void gsAffineOptParamMethod::reset(){
    m_obj = m_optParamMethod->evalObj();
    m_grad = m_optParamMethod->gradObj();
    m_hess = m_optParamMethod->hessObj(m_hessTagged);

    m_refFree = getFree();
    m_refTagged = getTagged();

    // For now we dont have any constraints so the KKT-system is simply the hessian
    m_KKTsystem = m_hess;

    // I problably have to inforce c_tagged = x, to allow change of tagge DoFs...

    m_rhs = -m_grad; // This rhs only works when the tagged has not been changed

    m_solver.compute(m_KKTsystem);
}

gsVector<> gsAffineOptParamMethod::getUpdate(gsVector<> x){
    // Setup rhs
    m_rhs = -m_grad - m_hessTagged*(x-m_refTagged); // This rhs only works when the tagged has not been changed

    // Solve optimality conditions to find an update to the free vars
    gsVector<> deltaFree = m_solver.solve(m_rhs);

    return m_refFree + deltaFree; // Return the new controlpoints ..

}

real_t gsAffineOptParamMethod::evalObj(gsVector<> c, gsVector<> x){
    return m_obj + m_grad.transpose()*(c-m_refFree) +
        0.5*(c-m_refFree).transpose()*m_hess*(c-m_refFree) +
        (x-m_refTagged).transpose()*m_hessTagged*(c-m_refFree);
}
