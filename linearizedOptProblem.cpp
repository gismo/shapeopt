#include <gismo.h>
#include <fstream>
#include "paraOptProblem.h"
#include "linearizedOptProblem.h"
using namespace gismo;

linearizedOptProblem::linearizedOptProblem(paraOptProblem* problem):m_problem(problem){
  m_cMatrix = problem->interfaceConstraintMatrix;
	index_t nfd = problem->iC.n_constraints;
  gsInfo << "nfd = " << nfd << "\n";
  m_cValues.setZero(nfd);

  obj = m_problem->evalObj();
  grad = m_problem->gradientObj();
  hess = m_problem->hessianObj();
  refCps = m_problem->getDesignVariables();

  nconst = m_cMatrix.rows();
  ndesign = m_cMatrix.cols();

  diffInBounds = m_problem->desUpperBounds() - m_problem->desLowerBounds();
  nz = (diffInBounds.array() == 0).count();

  M.setZero(nz,ndesign);

  index_t ind = 0;
  for(index_t i = 0; i < ndesign; i++){
      if (diffInBounds[i] == 0){
        M(ind,i) = 1;
        ind++;
      }

  }

  KKTsystem.setZero(nconst+ndesign+nz,nconst+ndesign+nz);

  KKTsystem.block(0,0,ndesign,ndesign) = hess;
  KKTsystem.block(0,ndesign,ndesign,nconst) = m_cMatrix.transpose();
  KKTsystem.block(ndesign,0,nconst,ndesign) = m_cMatrix;
  KKTsystem.block(ndesign+nconst,0,nz,ndesign) = M;
  KKTsystem.block(0,ndesign+nconst,ndesign,nz) = M.transpose();

  rhs.setZero(nconst + ndesign + nz);
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
