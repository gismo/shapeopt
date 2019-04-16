#include <gismo.h>
#include <fstream>
#include "harmonicOptProblem.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

harmonicOptProblem::harmonicOptProblem(gsMultiPatch<>* mpin):mp(mpin), dJC(mpin), iC(mpin){
  m_numDesignVars = dJC.n_controlpoints*2;
  m_numConstraints = iC.n_constraints;

  // Call getDvectors to generate LU factorization and stuff
  gsVector<> tmp(dJC.n_constraints);
  dJC.getDvectors(tmp);
  IpOptSparseMatrix J2 = iC.getJacobian();

  m_numConJacNonZero = J2.nnz();

  setDesignBounds(); //m_desLowerBounds and m_desUpperBounds is set in here.
  // gsInfo << "m_desUpperBounds  = \n" <<  m_desUpperBounds << "\n";

  // Set lower bounds to eps and 0's
  m_conUpperBounds.setZero(m_numConstraints);

  // gsInfo << "m_conLowerBounds  = \n" <<  m_conLowerBounds << "\n";

  // Set upper bounds to aBigNumber and 0's
  m_conLowerBounds.setZero(m_numConstraints);

  // gsInfo << "m_conUpperBounds  = " <<  m_conUpperBounds << "\n";

  // Concatenate J1 and J2 to get the final values
  m_conJacRows = J2.rows();
  m_conJacCols = J2.cols();

  m_curDesign = dJC.getDesignVariables();

  interfaceConstraintMatrix = iC.generateConstraintMatrix();
}

void harmonicOptProblem::setDesignBounds(){
  //FIXIT: Patch 3 is hardcoded in paraOptProblem
  if (mp->nBoxes() > 3){
    // gsInfo << "X vars on west boundaries freed ! (OBS: hardcoded)" << std::flush;
    // iC.freeBoundary(0,boundary::west,2);
    // iC.freeBoundary(1,boundary::west,2);
    // iC.freeBoundary(2,boundary::west,2);
    gsInfo << "X vars on south boundaries freed ! (OBS: hardcoded)";
    iC.freeBoundary(0,boundary::south,2);
    iC.freeBoundary(1,boundary::south,2);
    iC.freeBoundary(2,boundary::south,2);
    gsInfo << "Patch 3 is fixed (OBS: hardcoded)\n" << std::flush;
    m_desUpperBounds = iC.getUpperBounds(3); //Fix patch 3
  } else {
    gsInfo << "No patch is fixed (OBS: this happens when no. patches is less than 4)\n";
    m_desUpperBounds = iC.getUpperBounds(); //Fix no patch
  }
  m_desLowerBounds = iC.getLowerBounds(m_desUpperBounds);
}

gsVector<> harmonicOptProblem::getDesignVariables() const{
  return dJC.getDesignVariables();
}

void harmonicOptProblem::updateDesignVariables(gsVector<> des){
  return dJC.updateDesignVariables(des);
}

real_t harmonicOptProblem::evalObj() const{
  // gsInfo << "evalObj\n";
  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(*mp);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  // ev.options().setReal("quA",quA);
  // ev.options().setInt("quB",quB);
  // gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(*mp);

  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = ev.getVariable(x);
  variable fy = ev.getVariable(y);
  auto j00 = grad(fx)*jac(G)*grad(fx).tr();
  auto j10 = grad(fy)*jac(G)*grad(fx).tr();
  auto j01 = grad(fx)*jac(G)*grad(fy).tr();
  auto j11 = grad(fy)*jac(G)*grad(fy).tr();

  auto g11 = j00*j00 + j10*j10;
  auto g12 = j00*j01 + j10*j11;
  auto g22 = j01*j01 + j11*j11;

  gsFunctionExpr<> f1("x", "-y",2);
  variable ffun = ev.getVariable(f1);
  auto D = fjac(ffun);

  gsFunctionExpr<> fe("1.0","1.0",2);
  variable e = ev.getVariable(fe);

  return ev.integral((D*hess(G)*D*fform(G)).trace().sqNorm().val() + lambda_1 * hess(G).sqNorm() + lambda_2*jac(G).sqNorm());
  // return ev.integral((D*hess(G)*D*fform(G)).trace().sqNorm().val() + lambda_1 * (e.tr()*(hess(G)*hess(G)).trace()).val() + lambda_2*(jac(G)%jac(G)).val());// + lambda_2*jac(G)%jac(G));
  // return ev.integral(e.tr()*(hess(G)*hess(G)).trace());
}

real_t harmonicOptProblem::evalObj( const gsAsConstVector<real_t> & u) const {
  // gsInfo << "...evalObj\n";
  dJC.updateDesignVariables(u);
  return evalObj();
}

void harmonicOptProblem::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  result = interfaceConstraintMatrix*u;
}

void harmonicOptProblem::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  // gsInfo << "...jacobCon_into\n" << std::flush;
  dJC.updateDesignVariables(u);
  IpOptSparseMatrix J2 = iC.getJacobian();
  result = J2.values();
}

void harmonicOptProblem::print(){
  gsInfo << "m_numDesignVars  = " <<  m_numDesignVars << "\n";
  gsInfo << "m_numConstraints  = " <<  m_numConstraints << "\n";

  gsInfo << "m_numConJacNonZero  = " <<  m_numConJacNonZero << "\n";

  gsInfo << "m_conUpperBounds.size() = " << m_conUpperBounds.size() << "\n";
  gsInfo << "m_conLowerBounds.size() = " << m_conLowerBounds.size() << "\n";
  gsInfo << "m_desUpperBounds.size() = " << m_desUpperBounds.size() << "\n";
  gsInfo << "m_desLowerBounds.size() = " << m_desLowerBounds.size() << "\n";

  gsInfo << "m_curDesign.size() = " << m_curDesign.size() << "\n";

  // gsMatrix<> disp(m_desUpperBounds.size(),2);
  // disp << m_desUpperBounds,m_desLowerBounds;
  // gsInfo << ".. design upper and lower bounds\n";
  // gsInfo << disp << "\n";
  //
  gsMatrix<> disp2(m_conUpperBounds.size(),2);
  disp2 << m_conUpperBounds,m_conLowerBounds;
  gsInfo << ".. constraint upper and lower bounds\n";
  gsInfo << disp2 << "\n";

}

// void paraOptProblem::writeToFile(gsVector<> vec, std::string name) const{
//     gsInfo << "WRITING to " << name << "\n";
//   std::ofstream f(name);
//   for (auto &e : vec) f << std::setprecision(12) << e << "\n";
// }
//
// void paraOptProblem::writeToFile(gsMatrix<> mat, std::string name) const{
//     gsInfo << "WRITING to " << name << "\n";
//   std::ofstream f(name);
//   for(index_t i = 0; i < mat.rows(); i++){
//     for(index_t j = 0; j < mat.cols(); j++){
//       f << std::setprecision(20) << mat(i,j) << " ";
//     }
//     f << "\n";
//   }
// }
//
// void paraOptProblem::loadFromFile(std::string name){
//   std::ifstream f(name);
//   gsVector<> des(m_numDesignVars);
//   for(index_t i = 0; i < m_numDesignVars; i++){
//       f >> des[i];
//   }
//   dJC.updateDesignVariables(des);
//
// }
//
// gsMatrix<> paraOptProblem::hessianObj() const{
//   gsMatrix<> hess;
//   hess.setZero(m_numDesignVars,m_numDesignVars);
//   index_t ind = 0;
//   for(index_t i = 0; i < mp->nBoxes(); i++){
//     // gsInfo << "\ni = " << i << "\n";
//     gsMatrix<> xxMat,xyMat,yyMat;
//     evaluate2ndDerivOnPatch(i,xxMat,xyMat,yyMat);
//     index_t ncoefs = mp->patch(i).coefsSize();
//     hess.block(ind,ind,ncoefs,ncoefs) = xxMat;
//     hess.block(ind,ind+dJC.n_controlpoints,ncoefs,ncoefs) = xyMat.transpose();
//     hess.block(ind+dJC.n_controlpoints,ind,ncoefs,ncoefs) = xyMat;
//     hess.block(ind+dJC.n_controlpoints,ind+dJC.n_controlpoints,ncoefs,ncoefs) = yyMat;
//     ind += ncoefs;
//   }
//   return hess;
// }
//
// gsVector<> paraOptProblem::gradientObj() const{
//   gsVector<> grad;
//   grad.setZero(m_numDesignVars);
//   index_t ind = 0;
//   for(index_t i = 0; i < mp->nBoxes(); i++){
//     // gsInfo << "\ni = " << i << "\n";
//     gsVector<> xVec,yVec;
//     evaluateDerivOnPatch(i,xVec,yVec);
//     index_t ncoefs = mp->patch(i).coefsSize();
//     grad.segment(ind,ncoefs) = xVec;
//     grad.segment(ind+dJC.n_controlpoints,ncoefs) = yVec;
//     ind += ncoefs;
//   }
//   return grad;
// }

void harmonicOptProblem::reset(){
  setDesignBounds();
  m_curDesign = dJC.getDesignVariables();
}
