#include <gismo.h>
#include "maxDetJacOptProblem.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

maxDetJacOptProblem::maxDetJacOptProblem(gsMultiPatch<>* mpin):paraOptProblem(mpin), mp(mpin), dJC(mpin), iC(mpin){
  m_numDesignVars = dJC.n_controlpoints*mp->geoDim() + 1;
  m_numConstraints = dJC.n_constraints + iC.n_constraints;

  // Call getDvectors to generate LU factorization and stuff
  gsVector<> tmp(dJC.n_constraints);
  dJC.getDvectors(tmp);
  IpOptSparseMatrix J1 = -dJC.getJacobian();
  IpOptSparseMatrix J2 = iC.getJacobian();

  m_numConJacNonZero = J1.nnz() + J2.nnz() + dJC.n_constraints; // the constraint have a row of ones, the same length as the D vector

  setDesignBounds(); //m_desLowerBounds and m_desUpperBounds is set in here.
  // gsInfo << "m_desLowerBounds  = \n" <<  m_desLowerBounds.rows() << "\n";
  // gsInfo << "m_desUpperBounds  = \n" <<  m_desUpperBounds.rows() << "\n";

  // Set lower  bounds to 0's and upper bounds zero and a large number
  // we consider constraints -d_i > slack for all i
  m_conLowerBounds.setZero(m_numConstraints);

  m_conUpperBounds.setZero(m_numConstraints);
  m_conUpperBounds.segment(0,dJC.n_constraints).setOnes();
  m_conUpperBounds.segment(0,dJC.n_constraints) *= iC.aBigNumber;

  // gsInfo << "m_conLowerBounds  = \n" <<  m_conLowerBounds << "\n";
  // gsInfo << "m_conUpperBounds  = " <<  m_conUpperBounds << "\n";

  // Compute final part of gradient
  gsMatrix<> ones;
  ones.setZero(m_numConstraints,1);
  for (index_t i = 0; i < dJC.n_constraints; ++i){
    ones(i,0) = -1;
  }
  IpOptSparseMatrix J3(ones,0);


  // Concatenate J1 and J2
  J1.concatenate(J2,"col");

  // Concatenate J1 and J3 row wise to get the final values
  J1.concatenate(J3,"row");

  m_conJacRows = J1.rows();
  m_conJacCols = J1.cols();

  writeToFile(J1.asDense(),"../results/J1.txt");
  writeToFile(dJC.getJacobian().asDense(),"../results/dJ.txt");

  // Again implemented assuming negative determinant
  real_t minD = -dJC.getDvectors().maxCoeff();
  gsInfo << "\nmin D: " << minD << "\n";

  m_curDesign.setZero(m_numDesignVars,1);
  m_curDesign.block(0,0,m_numDesignVars-1,1) = dJC.getDesignVariables();
  m_curDesign(m_numDesignVars-1,0) = minD;

  interfaceConstraintMatrix = iC.generateConstraintMatrix();
}

void maxDetJacOptProblem::setDesignBounds(){
  m_desLowerBounds.setZero(m_numDesignVars);
  m_desUpperBounds.setZero(m_numDesignVars);

  //FIXIT: Patch 3 is (almost) hardcoded in paraOptProblem
  gsVector<> upperBounds;
  if (mp->nBoxes() > fixedPatch){
    gsInfo << "Patch 3 is fixed (OBS: hardcoded)\n" << std::flush;
    upperBounds = iC.getUpperBounds(fixedPatch); //Fix patch 3
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
    upperBounds  = iC.getUpperBounds(); //Fix no patch
  }
  m_desLowerBounds.segment(0,m_numDesignVars-1) = iC.getLowerBounds(upperBounds);
  m_desUpperBounds.segment(0,m_numDesignVars-1) = upperBounds;

  // No bound on slack variable (the last element)
  m_desLowerBounds[m_numDesignVars - 1] = -iC.aBigNumber;
  m_desUpperBounds[m_numDesignVars - 1] = iC.aBigNumber;
}

real_t maxDetJacOptProblem::evalObj ( const gsAsConstVector<real_t> & u) const {
  // gsInfo << "evalObj\n" << std::flush;
  return -u[m_numDesignVars - 1]; // Return slack variable
}

void maxDetJacOptProblem::gradObj_into(const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const{
  // gsInfo << "gradObj\n" << std::flush;
  // Set almost almost the most
  for (index_t i = 0; i < m_numDesignVars - 1; i++){
    result[i] = 0;
  }
  result[m_numDesignVars - 1] = -1;
}

void maxDetJacOptProblem::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  // gsInfo << "...evalCon_into" << "\n" << std::flush;

  gsVector<> des = u.segment(0,m_numDesignVars-1);
  dJC.updateDesignVariables(des);

  // we consider constraints -d_i > slack for all i
  result.segment(0,dJC.n_constraints) = -dJC.getDvectors();

  for (index_t i = 0; i < dJC.n_constraints; i++){
    result[i] -= u[m_numDesignVars - 1];
  }

  // gsInfo << "Min of constraint " << result.segment(0,dJC.n_constraints).minCoeff() << "\n";
  // gsInfo << "Max of constraint " << result.segment(0,dJC.n_constraints).maxCoeff() << "\n";


  result.segment(dJC.n_constraints,iC.n_constraints) = interfaceConstraintMatrix*des;

  // gsInfo << "max ifconst : " << result.segment(dJC.n_constraints,iC.n_constraints).norm() << "\n";
}

void maxDetJacOptProblem::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  // gsInfo << "...jacobCon_into\n" << std::flush;

  dJC.updateDesignVariables(u.segment(0,m_numDesignVars-1));

  // gsInfo << "min d: " << -dJC.getDvectors().maxCoeff() << "\n";
  // gsInfo << "slack variable : " << u[m_numDesignVars-1] << "\n";

  IpOptSparseMatrix J1 = -dJC.getJacobian();
  IpOptSparseMatrix J2 = iC.getJacobian();

  // Compute final part of gradient
  gsMatrix<> ones;
  ones.setZero(m_numConstraints,1);
  for (index_t i = 0; i < dJC.n_constraints; ++i){
    ones(i,0) = -1;
  }
  IpOptSparseMatrix J3(ones,0);

  // Concatenate J1 and J2
  J1.concatenate(J2,"col");

  // Concatenate J1 and J3 row wise to get the final values
  J1.concatenate(J3,"row");

  result = J1.values();
}

void maxDetJacOptProblem::reset(){
  setDesignBounds();

  // Again implemented assuming negative determinant
  real_t minD = -dJC.getDvectors().maxCoeff();

  m_curDesign.block(0,0,m_numDesignVars-1,1) = dJC.getDesignVariables();
  m_curDesign(m_numDesignVars-1,0) = minD;

  interfaceConstraintMatrix = iC.generateConstraintMatrix();
}

void maxDetJacOptProblem::writeToFile(gsVector<> vec, std::string name) const{
  std::ofstream f(name);
  for (auto &e : vec) f << std::setprecision(12) << e << "\n";
}

void maxDetJacOptProblem::writeToFile(gsMatrix<> mat, std::string name) const{
  std::ofstream f(name);
  for(index_t i = 0; i < mat.rows(); i++){
    for(index_t j = 0; j < mat.cols(); j++){
      f << mat(i,j) << " ";
    }
    f << "\n";
  }
}
