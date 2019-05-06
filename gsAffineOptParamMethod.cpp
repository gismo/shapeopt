#include <gismo.h>
#include "gsAffineParamMethod.h"
using namespace gismo;

gsAffineOptParamMethod::gsAffineOptParamMethod(gsOptParamMethod* pM): m_pM(pM)
{

  m_obj = m_optParamMethod->evalObj();
  m_grad = m_optParamMethod->gradObj();
  m_hess = m_optParamMethod->hessObj();

  m_refFree = getFree();
  m_refTagged = getTagged();

  // For now we dont have any constraints so the KKT-system is simply the hessian
  KKTsystem = hess;

  // I problably have to inforce c_tagged = x, to allow change of tagge DoFs...

  rhs.setZero();
  rhs.segment(0,ndesign) = -grad ;

  solver.compute(KKTsystem);
};

void linearizedOptProblem::reset(){
  obj = m_problem->evalObj();
  grad = m_problem->gradientObj();
  hess = m_problem->hessianObj();
  refCps = m_problem->getDesignVariables();

  // The only part of KKTsystem that depends on the reference parametrization are the hessian part
  KKTsystem.block(0,0,ndesign,ndesign) = hess;

  // Compute factorization of the new matrix
  solver.compute(KKTsystem);

  // Similarly rhs depends on grad
  rhs.segment(0,ndesign) = -grad ;

}

gsVector<> linearizedOptProblem::solve(gsVector<> deltaCps){

  // gsInfo << "Linearized obj BEFORE : " << evalObj(0,grad,hess,refCps,refCps) << "\n";
  // gsInfo << "Obj BEFORE : " << obj << "\n";

  gsVector<> Mrhs;
  Mrhs.setZero(nz);
  index_t ind = 0;
  for(index_t i = 0; i < ndesign; i++){
      if (diffInBounds[i] == 0){
        Mrhs(ind) = deltaCps[i];
        ind++;
      }

  }

  rhs.segment(ndesign + nconst,nz) = Mrhs;
  gsVector<> out = solver.solve(rhs);

  return out.segment(0,ndesign);
  // gsInfo << "Linearized obj AFTER : " << evalObj(0,grad,hess,m_problem->getDesignVariables(),refCps) << "\n";
  // gsInfo << "Liao obj AFTER : " << m_problem->evalObj() << "\n";

}

void linearizedOptProblem::solveAndUpdate(gsVector<> deltaCps){
  m_problem->updateDesignVariables(refCps + solve(deltaCps));
}

real_t linearizedOptProblem::evalObj(real_t obj, gsVector<> grad, gsMatrix<> hess, gsVector<> x, gsVector<> x0){
  return obj + grad.transpose()*(x-x0) + 0.5*(x-x0).transpose()*hess*(x-x0);
}
