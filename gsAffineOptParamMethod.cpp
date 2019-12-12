#include <gismo.h>
#include "gsAffineOptParamMethod.h"
using namespace gismo;

gsAffineOptParamMethod::gsAffineOptParamMethod(gsOptParamMethod* optParamMethod, bool use_Lagrangian):
    gsAffineParamMethod(optParamMethod->mp(),
    optParamMethod->mappers()),
    m_optParamMethod(optParamMethod),
    m_use_Lagrangian(use_Lagrangian)
{
    reset(); //Setup the problem with the parametrization hold on m_pM->m_mp as reference
};

void gsAffineOptParamMethod::reset(){
    if (m_use_Lagrangian){ // Use linearization of Lagrangian to define map
        gsInfo << "Use Lagrangian for linearization\n";
        m_obj = m_optParamMethod->evalLagrangian();
        m_grad = m_optParamMethod->gradLagrangian();
        m_hess = m_optParamMethod->hessLagrangian(m_hessTagged);
    } else { // Use linearization of objective (eg. Winslow) to define map
        gsInfo << "Use objective for linearization \n";
        m_obj = m_optParamMethod->evalObj();
        m_grad = m_optParamMethod->gradObj();
        m_hess = m_optParamMethod->hessObj(m_hessTagged);
    }
    gsInfo << "Norm of grad = " << m_grad.norm() << "\n\n";

    m_refFree = getFree();
    m_refTagged = getTagged();

    // For now we dont have any constraints so the KKT-system is simply the hessian
    m_KKTsystem = m_hess;

    m_rhs = -m_grad; // This rhs only works when the tagged has not been changed otherwise use: m_rhs = -m_grad - m_hessTagged*(x-m_refTagged);

    m_solver.compute(m_KKTsystem);
}

bool gsAffineOptParamMethod::updateAndReset()
{
    bool status = m_optParamMethod->update();
    reset();
    return status;
}

gsVector<> gsAffineOptParamMethod::getUpdate(gsVector<> x){
    // Setup rhs
    m_rhs = -m_grad - m_hessTagged*(x-m_refTagged);

    // Solve optimality conditions to find an update to the free vars
    gsVector<> deltaFree = m_solver.solve(m_rhs);

    return m_refFree + deltaFree; // Return the new controlpoints ..

}

real_t gsAffineOptParamMethod::evalObj(gsVector<> c, gsVector<> x){
    return m_obj + m_grad.transpose()*(c-m_refFree) +
        0.5*(c-m_refFree).transpose()*m_hess*(c-m_refFree) +
        (x-m_refTagged).transpose()*m_hessTagged*(c-m_refFree);
}
