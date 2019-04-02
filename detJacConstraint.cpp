#include <gismo.h>
#include "detJacConstraint.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

detJacConstraint::detJacConstraint(gsMultiPatch<>* mpin): mp(mpin), m_detJacBasis(*mpin){
  // Prepare basis for detJac
  // by setting the degree to (p-1)(p-1)
  m_isSolverSetup = false;
  int p = m_detJacBasis.maxCwiseDegree();
  m_detJacBasis.setDegree(2*p-1);
  // by reducing the continuity
  m_detJacBasis.reduceContinuity(1);
  gsInfo << "\n..detJacBasis.size: " << m_detJacBasis.size(0) << "\n";
  gsInfo << "\n..detJacBasis degree is: " << 2*p-1 << "\n";

  // Count the total number of controlpoints
  n_controlpoints = 0;
  for(int i = 0; i < mp->nBoxes(); i++){
    n_controlpoints += mp->patch(i).coefsSize();
    gsInfo << "Patch " << i << " has " << mp->patch(i).coefsSize() << " cps... \n";
  }

  n_constraints = m_detJacBasis.size();

  gsInfo << "n_controlpoints = " << n_controlpoints << "\n \n";
}

void detJacConstraint::getDvectors(gsVector<> &result){
  result = getDvectors();
}

gsVector<> detJacConstraint::getDvectors(){
  // gsInfo << "getDvectors\n" << std::flush;
  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  gsExprAssembler<> A(1,1);

  // Elements used for numerical integration
  A.setIntegrationElements(m_detJacBasis);
  gsExprEvaluator<> ev(A);

  space u = A.getSpace(m_detJacBasis);

  geometryMap G = A.getMap(*mp);

  A.initSystem();
  if (! m_isSolverSetup){
    A.assemble(u*u.tr());
    solverMassMatrix.compute(A.matrix());
    m_isSolverSetup = true;
  }

  gsMatrix<> solVector;
  solution u_sol = A.getSolution(u, solVector);

  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = A.getCoeff(x);
  variable fy = A.getCoeff(y);
  auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
  auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
  auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
  auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

  auto detJ = j00*j11 - j01*j10;

  A.assemble(u*detJ);

  solVector = solverMassMatrix.solve(A.rhs());

  // gsInfo << "solVector is a (" << solVector.rows() << ", " << solVector.cols() << ")\n";
  // gsInfo << "return!\n" << std::flush;
  return solVector;

}

IpOptSparseMatrix detJacConstraint::getJacobian(){
  // gsInfo << "getDetJacJacobian\n" << std::flush;
  // Prepare assembler
  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(*mp);
  A.setIntegrationElements(m_detJacBasis);
  gsExprEvaluator<> ev(A);

  // Define types
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(*mp);

  // Use simple functions to generate standard basis vectors [1,0] and [0,1]
  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = ev.getVariable(x);
  variable fy = ev.getVariable(y);

  // Setup and assemble the two matricies
  space u = A.getSpace(dbasis);
  space v = A.getTestSpace(u,m_detJacBasis);

  A.initMatrix();

  auto j00 = grad(fx)*jac(G)*grad(fx).tr();
  auto j10 = grad(fy)*jac(G)*grad(fx).tr();
  auto j01 = grad(fx)*jac(G)*grad(fy).tr();
  auto j11 = grad(fy)*jac(G)*grad(fy).tr();

  auto dxiR = grad(u)*grad(fx).tr();
  auto detaR = grad(u)*grad(fy).tr();

  auto v1 = dxiR*j11;
  auto v2 = -detaR*j10;

  auto z0gradu = v1 + v2;
  A.assemble(v*z0gradu.tr());

  gsMatrix<> xRhs = A.matrix();

  A.initSystem();

  auto v3 = -dxiR*j01;
  auto v4 = detaR*j00;

  auto z1gradu = v3 + v4;
  A.assemble(v*z1gradu.tr());

  gsMatrix<> yRhs = A.matrix();

  gsMatrix<> xJac = solverMassMatrix.solve(xRhs);
  gsMatrix<> yJac = solverMassMatrix.solve(yRhs);

  // FIXIT : How to exploit sparsity!!!??? TODO
  IpOptSparseMatrix xJ(xJac,-1);
  IpOptSparseMatrix yJ(yJac,-1);

  xJ.concatenate(yJ,"row");       // Store full jacobian in xMat
  return xJ;

}

gsVector<> detJacConstraint::getDesignVariables(){
  gsVector<> cx(n_controlpoints);
  gsVector<> cy(n_controlpoints);
  index_t j = 0;
  for(index_t i = 0; i < mp->nBoxes(); i++){
    for(index_t k = 0; k < mp->patch(i).coefsSize(); k++){
      cx[j] = mp->patch(i).coef(k,0);
      cy[j] = mp->patch(i).coef(k,1);
      j++;
    }
  }

  gsVector<> out(2*n_controlpoints);
  out << cx,cy;

  return out;
}

void detJacConstraint::updateDesignVariables(gsVector<> des){

  GISMO_ASSERT(des.size() == 2*n_controlpoints, "Design vector is of wrong size.!");

  gsVector<> cx = des.segment(0,n_controlpoints);
  gsVector<> cy = des.segment(n_controlpoints,n_controlpoints);

  index_t seg = 0;
  for(index_t i = 0; i < mp->nBoxes(); i++){
    index_t n_coefs = mp->patch(i).coefsSize();
    gsMatrix<> cc(n_coefs,2);

    cc << cx.segment(seg,n_coefs),
              cy.segment(seg,n_coefs);

    seg += n_coefs;

    mp->patch(i).setCoefs(cc);
  }
}

//FIXIT: detJacConstraint is coded for negative determinant.. FIXIT
gsVector<> detJacConstraint::getUpperBounds(real_t eps){
    gsVector<> out;
    out.setOnes(n_constraints);
    out *= -eps;
    return out;
  }

void detJacConstraint::plotDetJ(std::string name){
  // gsInfo << "getDvectors\n" << std::flush;
  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  gsExprAssembler<> A(1,1);

  // Elements used for numerical integration
  A.setIntegrationElements(m_detJacBasis);
  gsExprEvaluator<> ev(A);

  space u = A.getSpace(m_detJacBasis);

  geometryMap G = A.getMap(*mp);

  A.initSystem();
  if (! m_isSolverSetup){
    A.assemble(u*u.tr());
    solverMassMatrix.compute(A.matrix());
    m_isSolverSetup = true;
  }

  gsMatrix<> solVector;
  solution u_sol = A.getSolution(u, solVector);

  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = A.getCoeff(x);
  variable fy = A.getCoeff(y);
  auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
  auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
  auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
  auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

  auto detJ = j00*j11 - j01*j10;

  A.assemble(u*detJ);

  solVector = solverMassMatrix.solve(A.rhs());

  gsMultiPatch<> dJ;
  u_sol.extract(dJ);

  variable out = A.getCoeff(dJ);

	gsInfo<<"Plotting " << name << " in Paraview...\n";
  ev.writeParaview( out   , G, name);
	// ev.options().setSwitch("plot.elements", true);

}
