#include <gismo.h>
#include <fstream>
#include "paraOptProblem.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

paraOptProblem::paraOptProblem(gsMultiPatch<>* mpin):mp(mpin), dJC(mpin), iC(mpin){
  m_numDesignVars = dJC.n_controlpoints*2;
  m_numConstraints = dJC.n_constraints + iC.n_constraints;

  // Call getDvectors to generate LU factorization and stuff
  gsVector<> tmp(dJC.n_constraints);
  dJC.getDvectors(tmp);
  IpOptSparseMatrix J1 = dJC.getJacobian();
  IpOptSparseMatrix J2 = iC.getJacobian();

  m_numConJacNonZero = J1.nnz() + J2.nnz();

  setDesignBounds(); //m_desLowerBounds and m_desUpperBounds is set in here.
  // gsInfo << "m_desUpperBounds  = \n" <<  m_desUpperBounds << "\n";

  // Set lower bounds to eps and 0's
  m_conUpperBounds.setZero(m_numConstraints);
  // gsInfo << "getLowerBounds.size = " << dJC.getLowerBounds(0.1).size() << "\n";
  // gsInfo << "dJC.n_constraints = " << dJC.n_constraints << "\n";
  m_conUpperBounds.segment(0,dJC.n_constraints) = dJC.getUpperBounds(m_eps);

  // gsInfo << "m_conLowerBounds  = \n" <<  m_conLowerBounds << "\n";

  // Set upper bounds to aBigNumber and 0's
  m_conLowerBounds.setZero(m_numConstraints);
  m_conLowerBounds.segment(0,dJC.n_constraints).setOnes();
  m_conLowerBounds.segment(0,dJC.n_constraints) *= -iC.aBigNumber;

  // gsInfo << "m_conUpperBounds  = " <<  m_conUpperBounds << "\n";

  // Concatenate J1 and J2 to get the final values
  J1.concatenate(J2,"col");
  m_conJacRows = J1.rows();
  m_conJacCols = J1.cols();

  m_curDesign = dJC.getDesignVariables();

  interfaceConstraintMatrix = iC.generateConstraintMatrix();
}

void paraOptProblem::setDesignBounds(){
  //FIXIT: Patch 3 is hardcoded in paraOptProblem
  if (mp->nBoxes() > 3){
    gsInfo << "Patch 3 is fixed (OBS: hardcoded)\n" << std::flush;
    m_desUpperBounds = iC.getUpperBounds(3); //Fix patch 3
    // gsInfo << "X vars on west boundaries freed ! (OBS: hardcoded)" << std::flush;
    // iC.freeBoundary(0,boundary::west,2);
    // iC.freeBoundary(1,boundary::west,2);
    // iC.freeBoundary(2,boundary::west,2);
    gsInfo << "X vars on south boundaries freed ! (OBS: hardcoded)";
    iC.freeBoundary(0,boundary::south,2);
    iC.freeBoundary(1,boundary::south,2);
    iC.freeBoundary(2,boundary::south,2);
  } else {
    gsInfo << "No patch is fixed (OBS: this happens when no. patches is less than 4)\n";
    m_desUpperBounds = iC.getUpperBounds(); //Fix no patch
  }
  m_desLowerBounds = iC.getLowerBounds(m_desUpperBounds);
}

gsVector<> paraOptProblem::getDesignVariables() const{
  return dJC.getDesignVariables();
}

void paraOptProblem::updateDesignVariables(gsVector<> des){
  return dJC.updateDesignVariables(des);
}

real_t paraOptProblem::evalObj() const{
  // gsInfo << "evalObj\n";
  real_t result = 0;
  for(index_t i = 0; i < mp->nBoxes(); i++){
    result += evaluateOnPatch(i);
  }
  return result;
}

void paraOptProblem::gradObj_into(const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const{
  // gsInfo << "gradObj\n";
  // gsInfo << "\n..gradObj_into\n";
  dJC.updateDesignVariables(u);
  index_t ind = 0;
  for(index_t i = 0; i < mp->nBoxes(); i++){
    // gsInfo << "\ni = " << i << "\n";
    gsVector<> xVec,yVec;
    evaluateDerivOnPatch(i,xVec,yVec);
    index_t ncoefs = mp->patch(i).coefsSize();
    result.segment(ind,ncoefs) = xVec;
    result.segment(ind+dJC.n_controlpoints,ncoefs) = yVec;
    ind += ncoefs;
  }
}

real_t paraOptProblem::evalObj( const gsAsConstVector<real_t> & u) const {
  // gsInfo << "...evalObj\n";
  dJC.updateDesignVariables(u);
  return evalObj();
}

void paraOptProblem::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  // gsInfo << "...evalCon_into" << "\n";
  gsVector<> tmp(dJC.n_constraints);

  // gsInfo << "n_constraints : " << dJC.n_constraints << "\n";

  dJC.updateDesignVariables(u);
  dJC.getDvectors(tmp);

  // gsInfo << "max d: " << tmp.maxCoeff() << "\n";

  result.segment(0,dJC.n_constraints) = tmp;
  result.segment(dJC.n_constraints,iC.n_constraints) = interfaceConstraintMatrix*u;

  // gsInfo << "max ifconst : " << result.segment(dJC.n_constraints,iC.n_constraints).norm() << "\n";
}

void paraOptProblem::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  // gsInfo << "...jacobCon_into\n" << std::flush;
  dJC.updateDesignVariables(u);
  IpOptSparseMatrix J1 = dJC.getJacobian();
  IpOptSparseMatrix J2 = iC.getJacobian();
  J1.concatenate(J2,"col");
  result = J1.values();
}

void paraOptProblem::print(){
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

void paraOptProblem::writeToFile(gsVector<> vec, std::string name) const{
  std::ofstream f(name);
  for (auto &e : vec) f << std::setprecision(12) << e << "\n";
}

void paraOptProblem::writeToFile(gsMatrix<> mat, std::string name) const{
  std::ofstream f(name);
  for(index_t i = 0; i < mat.rows(); i++){
    for(index_t j = 0; j < mat.cols(); j++){
      f << mat(i,j) << " ";
    }
    f << "\n";
  }
}

void paraOptProblem::loadFromFile(std::string name){
  std::ifstream f(name);
  gsVector<> des(m_numDesignVars);
  for(index_t i = 0; i < m_numDesignVars; i++){
      f >> des[i];
  }
  dJC.updateDesignVariables(des);

}

gsMatrix<> paraOptProblem::hessianObj() const{
  gsMatrix<> hess;
  hess.setZero(m_numDesignVars,m_numDesignVars);
  index_t ind = 0;
  for(index_t i = 0; i < mp->nBoxes(); i++){
    // gsInfo << "\ni = " << i << "\n";
    gsMatrix<> xxMat,xyMat,yyMat;
    evaluate2ndDerivOnPatch(i,xxMat,xyMat,yyMat);
    index_t ncoefs = mp->patch(i).coefsSize();
    hess.block(ind,ind,ncoefs,ncoefs) = xxMat;
    hess.block(ind,ind+dJC.n_controlpoints,ncoefs,ncoefs) = xyMat.transpose();
    hess.block(ind+dJC.n_controlpoints,ind,ncoefs,ncoefs) = xyMat;
    hess.block(ind+dJC.n_controlpoints,ind+dJC.n_controlpoints,ncoefs,ncoefs) = yyMat;
    ind += ncoefs;
  }
  return hess;
}

gsVector<> paraOptProblem::gradientObj() const{
  gsVector<> grad;
  grad.setZero(m_numDesignVars);
  index_t ind = 0;
  for(index_t i = 0; i < mp->nBoxes(); i++){
    // gsInfo << "\ni = " << i << "\n";
    gsVector<> xVec,yVec;
    evaluateDerivOnPatch(i,xVec,yVec);
    index_t ncoefs = mp->patch(i).coefsSize();
    grad.segment(ind,ncoefs) = xVec;
    grad.segment(ind+dJC.n_controlpoints,ncoefs) = yVec;
    ind += ncoefs;
  }
  return grad;
}

void paraOptProblem::reset(){
  setDesignBounds();
  m_curDesign = dJC.getDesignVariables();
}
