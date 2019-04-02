/** @file main.cpp

    @brief Test of calculating derivatives of detJ

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Asger Limkilde
*/

#include <gismo.h>
#include "detJacConstraint.h"

using namespace gismo;

real_t evalObj(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
  dJC.updateDesignVariables(des);

  gsExprAssembler<> A(1,1);
  gsExprEvaluator<> ev(A);

  A.options().setInt("quB",quB);
  A.options().setReal("quA",quA);
  ev.options().setInt("quB",quB);
  ev.options().setReal("quA",quA);

  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  gsMultiBasis<> dbasis(mp);
  A.setIntegrationElements(dbasis);

  geometryMap G = A.getMap(mp);

  return ev.integral(jac(G).det());
}

gsVector<> gradObjImpl1(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
  dJC.updateDesignVariables(des);

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(mp);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  A.options().setInt("quB",quB);
  A.options().setReal("quA",quA);
  ev.options().setInt("quB",quB);
  ev.options().setReal("quA",quA);

  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(mp);

  space u = A.getSpace(dbasis);


  A.initSystem();
  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = A.getCoeff(x);
  variable fy = A.getCoeff(y);
  auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
  auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
  auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
  auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

  auto detJ = j00*j11 - j10*j01;
	auto detJinv = 1/detJ.val();

	auto signOfDetJ = detJinv*meas(G);

  auto d_detJ_dcx = signOfDetJ*(uxi*j11 - ueta*j10) ;
  auto d_detJ_dcy = signOfDetJ*(ueta*j00 - uxi*j01) ;

  A.assemble(d_detJ_dcx);

  gsMatrix<> xVec = A.rhs();

  A.initSystem();

  A.assemble(d_detJ_dcy);

  gsMatrix<> yVec = A.rhs();

  gsMatrix<> out(xVec.rows()*2,1);
  out << xVec, yVec;

  return out;
}

gsVector<> gradObjImpl2(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
  dJC.updateDesignVariables(des);

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(mp);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  A.options().setInt("quB",quB);
  A.options().setReal("quA",quA);
  ev.options().setInt("quB",quB);
  ev.options().setReal("quA",quA);

  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(mp);

  space u = A.getSpace(dbasis);


  A.initSystem();
  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = A.getCoeff(x);
  variable fy = A.getCoeff(y);
  auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
  auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
  auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
  auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

  auto detJ = j00*j11 - j10*j01;
	auto detJinv = 1/detJ.val();

	auto signOfDetJ = detJ.sgn();

  auto d_detJ_dcx = signOfDetJ*(uxi*j11 - ueta*j10) ;
  auto d_detJ_dcy = signOfDetJ*(ueta*j00 - uxi*j01) ;

  A.assemble(d_detJ_dcx);

  gsMatrix<> xVec = A.rhs();

  A.initSystem();

  A.assemble(d_detJ_dcy);

  gsMatrix<> yVec = A.rhs();

  gsMatrix<> out(xVec.rows()*2,1);
  out << xVec, yVec;

  return out;
}

gsVector<> gradObjImpl3(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
  dJC.updateDesignVariables(des);

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(mp);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  A.options().setInt("quB",quB);
  A.options().setReal("quA",quA);
  ev.options().setInt("quB",quB);
  ev.options().setReal("quA",quA);

  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(mp);

  space u = A.getSpace(dbasis,mp.geoDim()); // We are trying to find a vector expression

  A.initSystem();
  A.assemble(jac(G).det()*matrix_by_space(jac(G).inv(),jac(u)).trace());

  return A.rhs();
}

void convergenceTestOfJacobian(real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC,index_t impl){
	gsVector<> des = dJC.getDesignVariables();
  dJC.updateDesignVariables(des);
  index_t nd = des.size();

	// std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(nd);

	gsInfo << "\n Size of design vector : " << nd << "\n";

	real_t obj = evalObj(des,quA,quB,mp,dJC);
  gsVector<> grad;
  if (impl == 1){
	   grad = gradObjImpl1(des,quA,quB,mp,dJC);
   } else if (impl == 2) {
	   grad = gradObjImpl2(des,quA,quB,mp,dJC);
   } else if (impl == 3) {
	   grad = gradObjImpl2(des,quA,quB,mp,dJC);
   } else {
     GISMO_ERROR("Wrong impl choice!!!\n");
   }

	index_t beg = 4;
	index_t n = 20;
	gsVector<> Eps(n);
	gsVector<> Error0(n);
	gsVector<> Error1(n);

	for(index_t i = 0; i < n; i++){
		// Generate pertubation
		real_t eps = pow(2,-beg-i);
		Eps[i] = eps;

		gsVector<> perturp = ran;

		perturp /= perturp.norm();

		perturp *= eps;

		gsVector<> newDes = des + perturp;

		real_t newObj = evalObj(newDes,quA,quB,mp,dJC);
  	real_t guess0 = obj;
		real_t guess1 = obj + grad.transpose()*perturp;

		gsInfo << guess0 <<" " << guess1 << " " << newObj << "\n";

		real_t error0 = std::abs(guess0 - newObj);
		real_t error1 = std::abs(guess1 - newObj);

		Error0[i] = error0;
		Error1[i] = error1;
	}

	gsVector<> rate;
	rate.setZero(n);
	rate.segment(1,n-1) = log10(Error1.segment(1,n-1).array()/Error1.segment(0,n-1).array())/log10(2);
	gsMatrix<> disp(n,4);
	disp << Eps,Error0,Error1,rate;
	gsInfo << "\neps \tErr0 \tErr1 \trate";
	gsInfo << disp << "\n";

}

int main(int argc, char* argv[])
{

  // Parse command line
  std::string output("");
  int degree = 2;
  int numRefine = 1;
  bool plot = false;

  int quA = 1;
  int quB = 1;

  int impl = 1;

  gsCmdLine cmd("A test of lumped mass matricies");
  cmd.addInt("p", "degree", "Degree of B-Splines.", degree);
  cmd.addInt("n", "numberRefine", "Number of refinements", numRefine);
  cmd.addInt("i", "howToCalculateDeriv", "How to calculate derivatives", impl);

  cmd.addInt("A", "quA", "quA", quA);
  cmd.addInt("B", "quB", "quB", quB);

  cmd.addString("o","output","Name of the output file",output);
  cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

  cmd.getValues(argc,argv);

  char buffer [50];
  std::sprintf(buffer,"p = %d \n n = %d\n",degree,numRefine);
  gsInfo << buffer;

  gsMultiPatch<> mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
  mp.basis(0).setDegree(degree);

  for(int i = 0; i < numRefine; i++){
    mp.uniformRefine();
  }

  int len = mp[0].coefsSize();

  // Change direction to test sign
  // int n = sqrt(len);
  // gsInfo << "n = " << n << " and len = " << len << "\n";
  //
  // gsMatrix<> cc = mp[0].coefs();
  // gsMatrix<> cx = cc.block(0,0,len,1);
  // gsMatrix<> cy = cc.block(0,1,len,1);
  //
  // cx.resize(n,n);
  // cy.resize(n,n);
  //
  // for (index_t i = 0; i < n; i++){
  //   cx.col(i).reverseInPlace();
  //   cy.col(i).reverseInPlace();
  // }
  //
  // cx.resize(len,1);
  // cy.resize(len,1);
  //
  // cc.block(0,0,len,1) = cx;
  // cc.block(0,1,len,1) = cy;
  //
  // mp.patch(0).setCoefs(cc);

  detJacConstraint dJC(&mp);
  gsVector<> des = dJC.getDesignVariables();

  // gsInfo << dJC.getDvectors();
  // real_t eps = 0.5;
  //
	// std::srand((unsigned int) std::time(0));
  // gsVector<> ran;
	// ran.setRandom(des.size());
  // ran *= eps;

  gsInfo << "The domain is a "<< mp <<"\n";

  gsInfo << "evalObj() = " << evalObj(des,quA,quB,mp,dJC);
  gsVector<> grad1 = gradObjImpl1(des,quA,quB,mp,dJC);
  gsVector<> grad2 = gradObjImpl2(des,quA,quB,mp,dJC);
  gsVector<> grad3 = gradObjImpl3(des,quA,quB,mp,dJC);

  gsInfo << "Difference in implementation 1 and 2 : " << (grad1-grad2).norm() << "\n";

  gsMatrix<> disp(grad1.rows(),3);
  disp << grad1, grad2, grad3;
  gsInfo << "\n" << disp << "\n" ;

  // convergenceTestOfJacobian(quA,quB,mp,dJC,impl);

  if (!output.empty()){
    gsInfo << "Output to " << output << "\n";
  }

  gsInfo << "\n\n==== DONE ====\n";
  return 0;
}
