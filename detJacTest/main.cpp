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

  auto d_detJ_dcx = signOfDetJ.val()*(uxi*j11.val() - ueta*j10.val()) ;
  auto d_detJ_dcy = signOfDetJ.val()*(ueta*j00.val() - uxi*j01.val()) ;

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

  auto detJ = j00.val()*j11.val() - j10.val()*j01.val();
	auto detJinv = 1/detJ.val();

	auto signOfDetJ = detJ.sgn();

  auto d_detJ_dcx = signOfDetJ*(uxi*j11.val() - ueta*j10.val()) ;
  auto d_detJ_dcy = signOfDetJ*(ueta*j00.val() - uxi*j01.val()) ;

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
  A.assemble(jac(G).det()*
  matrix_by_space(jac(G).inv(),jac(u)).trace()
  // (jac(u)*jac(G).inv()).trace()
  // (- jac(G).inv() *..under construction.. jac(u) * jac(G).inv())
  );

  return A.rhs();
}

real_t evalObj2(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
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

  gsFunctionExpr<> f("x + 2*y",2);
  variable ff = A.getCoeff(f);

  geometryMap G = A.getMap(mp);

  return ev.integral((fjac(ff).tr()*jac(G).inv())*jac(G).inv().tr()*fjac(ff));
}

gsVector<> gradObjImpl12(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
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

  gsFunctionExpr<> f("x + 2*y",2);
  variable ff = A.getCoeff(f);

  A.initSystem();
  A.assemble(-collapse(matrix_by_space_tr(jac(G).inv(),jac(u)),(jac(G).inv().tr()*fjac(ff))) * fjac(ff));
  // (jac(u)*jac(G).inv()).trace()
  // (- jac(G).inv() *..under construction.. jac(u) * jac(G).inv())

  return A.rhs();
}

gsVector<> gradObjImpl22(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
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

  gsFunctionExpr<> f("x + 2*y",2);
  variable ff = A.getCoeff(f);

  gsInfo << "test\n";

  A.initSystem();
  A.assemble(-collapse(fjac(ff).tr(),matrix_by_space(jac(G).inv(),jac(u))*jac(G).inv()) * jac(G).inv().tr()*fjac(ff)
    + -collapse(fjac(ff).tr()*jac(G).inv(),matrix_by_space_tr(jac(G).inv(),jac(u)))*jac(G).inv().tr() * fjac(ff)
    );
  // (jac(u)*jac(G).inv()).trace()
  // (- jac(G).inv() *..under construction.. jac(u) * jac(G).inv())
  return A.rhs();
}

gsVector<> evalVec(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
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

  space u = A.getSpace(dbasis);

  gsFunctionExpr<> f("x + 2*y",2);
  variable ff = A.getCoeff(f);

  geometryMap G = A.getMap(mp);

  A.initSystem();
  // A.assemble(igrad(u,G)*fjac(ff)*jac(G).det());
  // A.assemble(igrad(u,G)*(jac(G).inv().tr()*fjac(ff))*jac(G).det());
  // A.assemble(igrad(u,G)*jac(G).inv().tr()*fjac(ff)*jac(G).det());
  A.assemble(igrad(u,G)*jac(G).inv().tr()*fjac(ff)*meas(G));

  return A.rhs();
}

gsMatrix<> gradVec(gsVector<> des,real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC){
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

  space v = A.getTestSpace(u,dbasis,mp.geoDim()); // We are trying to find a vector expression

  gsFunctionExpr<> f("x + 2*y",2);
  variable ff = A.getCoeff(f);

  gsInfo << "test\n";

  A.initSystem();
  // A.assemble(matrix_by_space(jac(G).inv(),jac(v)).trace()*fjac(ff).tr()*igrad(u,G).tr()*jac(G).det()
    // - collapse(fjac(ff).tr(),matrix_by_space(jac(G).inv(),jac(v))*jac(G).inv()) * igrad(u,G).tr()*jac(G).det()
    // - collapse(fjac(ff).tr()*jac(G).inv(),matrix_by_space_tr(jac(G).inv(),jac(v)))*igrad(u,G).tr() * jac(G).det()
    // );
  // A.assemble(matrix_by_space(jac(G).inv(),jac(v)).trace()*(igrad(u,G)*jac(G).inv().tr()*fjac(ff)).tr()*meas(G)
  A.assemble(matrix_by_space(jac(G).inv(),jac(v)).trace()*fjac(ff).tr()*jac(G).inv()*igrad(u,G).tr()*meas(G)
    - collapse(fjac(ff).tr(),matrix_by_space(jac(G).inv(),jac(v))*jac(G).inv()) * igrad(u,G).tr()*meas(G)
    - collapse(fjac(ff).tr()*jac(G).inv(),matrix_by_space_tr(jac(G).inv(),jac(v)))*igrad(u,G).tr() *meas(G)
  // A.assemble(matrix_by_space(jac(G).inv(),jac(v)).trace()*(grad(u)*jac(G).inv().tr()*fjac(ff)).tr()*jac(G).det()
  //   - collapse(fjac(ff).tr(),matrix_by_space(jac(G).inv(),jac(v))*jac(G).inv()) * igrad(u,G).tr()*jac(G).det()
  //   - collapse(fjac(ff).tr()*jac(G).inv(),matrix_by_space(jac(G).inv(),jac(v))*jac(G).inv())*grad(u).tr() * jac(G).det()
);
  // (jac(u)*jac(G).inv()).trace()
  // (- jac(G).inv() *..under construction.. jac(u) * jac(G).inv())
  return A.matrix();

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
	   grad = gradObjImpl3(des,quA,quB,mp,dJC);
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

void convergenceTestOfJacobian2(real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC,index_t impl){
    gsVector<> des = dJC.getDesignVariables();
  dJC.updateDesignVariables(des);
  index_t nd = des.size();

	// std::srand((unsigned int) std::time(0));
	gsVector<> ran;
	ran.setRandom(nd);

	gsInfo << "\n Size of design vector : " << nd << "\n";

	real_t obj = evalObj2(des,quA,quB,mp,dJC);
    gsVector<> grad;
    if (impl == 1)
        grad = gradObjImpl12(des,quA,quB,mp,dJC);
    if (impl == 2)
        grad = gradObjImpl22(des,quA,quB,mp,dJC);

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

		real_t newObj = evalObj2(newDes,quA,quB,mp,dJC);
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

void convergenceTestOfJacobian3(real_t quA,int quB,gsMultiPatch<> &mp,detJacConstraint &dJC,index_t impl){
    gsVector<> des = dJC.getDesignVariables();
    dJC.updateDesignVariables(des);
    index_t nd = des.size();

    // std::srand((unsigned int) std::time(0));
    gsVector<> ran;
    ran.setRandom(nd);

	gsInfo << "\n Size of design vector : " << nd << "\n";

	gsVector<> vec = evalVec(des,quA,quB,mp,dJC);
    gsMatrix<> grad = gradVec(des,quA,quB,mp,dJC);
    gsInfo << "grad: (" << grad.rows() << ", " << grad.cols() << ")\n";

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

		gsVector<> newVec = evalVec(newDes,quA,quB,mp,dJC);
  	    gsVector<> guess0 = vec;
		gsVector<> guess1 = vec + grad.transpose()*perturp;

		// gsInfo << guess0 <<" " << guess1 << " " << newObj << "\n";

		real_t error0 = (guess0 - newVec).norm();
		real_t error1 = (guess1 - newVec).norm();

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

void readFromTxt(std::string name, gsMatrix<> &matrix){
	std::ifstream infile;
	infile.open(name);
	gsInfo << "Loading from " << name << "\n";
	for(int i = 0; i < matrix.rows(); i++){
		for(int j = 0; j < matrix.cols(); j++){
			infile >> matrix(i,j);
		}
	}
	infile.close();
}

gsVector<> loadVec(index_t n,std::string name){
		gsVector<> vec;
		vec.setZero(n);
		std::ifstream file (name);
		real_t val;
		for(index_t i = 0; i < n; i++){
			file >> val;
			vec[i] = val;
		}
		return vec;

}

gsMultiPatch<> getGeometry(index_t n, index_t m, index_t degree){

    std::string folder;
    if (n == 4 && m == 4 && degree == 2){
        folder = "/../parametrizations/para_x4_y4_p2_q2/";
    } else if (n == 5 && m == 4 && degree == 2){
        folder = "/../parametrizations/para_x5_y4_p2_q2/";
    } else if (n == 6 && m == 4 && degree == 2){
        folder = "/../parametrizations/para_x6_y4_p2_q2/";
    } else if (n == 7 && m == 4 && degree == 2){
        folder = "/../parametrizations/para_x7_y4_p2_q2/";
    } else if (n == 8 && m == 4 && degree == 2){
        folder = "/../parametrizations/para_x8_y4_p2_q2/";
    } else if (n == 4 && m == 5 && degree == 2){
        folder = "/../parametrizations/para_x4_y5_p2_q2/";
    } else if (n == 5 && m == 5 && degree == 2){
        folder = "/../parametrizations/para_x5_y5_p2_q2/";
    } else if (n == 6 && m == 5 && degree == 2){
        folder = "/../parametrizations/para_x6_y5_p2_q2/";
    } else if (n == 7 && m == 5 && degree == 2){
        folder = "/../parametrizations/para_x7_y5_p2_q2/";
    } else if (n == 8 && m == 5 && degree == 2){
        folder = "/../parametrizations/para_x8_y5_p2_q2/";
    }  else if (n == 6 && m == 6 && degree == 2){
        folder = "/../parametrizations/para_x6_y6_p2_q2/";
    } else {
        GISMO_ERROR("Parametrization is not generated with these parameters..\n");
    }



	gsInfo << "----------------------\n\n"
	<< "n: " << n << "\n\n"
	<< "m: " << m << "\n\n"
	<< "degree: " << degree << "\n\n"
	<< "----------------------\n\n";

	// 1. construction of a knot vector for each direction
	gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);
	// 2. construction of a basis
	gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
	// 3. construction of a coefficients
	gsMatrix<> greville = basis.anchors();
	gsMatrix<> coefs (greville.cols(), 2);

	readFromTxt(BASE_FOLDER + folder + "l.txt", coefs);

	gsInfo << coefs;

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  left(basis, coefs);

	//! [Geometry data]
	// Define Geometry, must be a gsMultiPatch object

	gsMultiPatch<> patches = gsMultiPatch<>(left);

	readFromTxt(BASE_FOLDER + folder + "b.txt", coefs);

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  bottom(basis, coefs);
	patches.addPatch(bottom);

	readFromTxt(BASE_FOLDER + folder + "r.txt", coefs);

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  right(basis, coefs);
	patches.addPatch(right);

	gsTensorBSplineBasis<2, real_t> basisMid(kv1, kv1);
	// 3. construction of a coefficients
	gsMatrix<> grevilleMid = basisMid.anchors();
	gsMatrix<> coefsMid (grevilleMid.cols(), 2);

	readFromTxt(BASE_FOLDER + folder + "m.txt", coefsMid);
	gsTensorBSpline<2, real_t>  middle(basisMid, coefsMid);
	patches.addPatch(middle);

	readFromTxt(BASE_FOLDER + folder + "t.txt", coefs);

	// 4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  top(basis, coefs);
	patches.addPatch(top);

	double tol = 1e-2;
	patches.computeTopology(tol,true);
	patches.closeGaps(tol);

	std::string out = "Geometry";
	gsInfo << "Writing the gsMultiPatch to a paraview file: " << out
	<< "\n\n";
	gsWriteParaview(patches, out);

	return patches;
	// GISMO_ERROR("stop...");

}

int main(int argc, char* argv[])
{

  // Parse command line
  std::string output("");
  int degree = 2;
  int numRefine = 1;
  bool plot = false;
  int startDes = -1;
  int quA = 1;
  int quB = 1;

  int impl = 1;

  gsCmdLine cmd("A test of lumped mass matricies");
  cmd.addInt("p", "degree", "Degree of B-Splines.", degree);
  cmd.addInt("n", "numberRefine", "Number of refinements", numRefine);
  cmd.addInt("i", "howToCalculateDeriv", "How to calculate derivatives", impl);
  cmd.addInt("s", "startDes", "startDes", startDes);

  cmd.addInt("A", "quA", "quA", quA);
  cmd.addInt("B", "quB", "quB", quB);

  cmd.addString("o","output","Name of the output file",output);
  cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

  cmd.getValues(argc,argv);

  char buffer [50];
  std::sprintf(buffer,"p = %d \n n = %d\n",degree,numRefine);
  gsInfo << buffer;

  // gsMultiPatch<> mp = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));
  // mp.basis(0).setDegree(degree);
  // gsMultiPatch<> mp(getGeometry(5,4,2).patch(1));
  gsMultiPatch<> mp = getGeometry(5,4,2);

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

  std::srand((unsigned int) std::time(0));
  detJacConstraint dJC(&mp);
  gsVector<> des = dJC.getDesignVariables();
  gsVector<> perturp;
  perturp.setRandom(des.size());

  perturp /= perturp.norm();
  perturp *= 0.4;

  // des += perturp;
  // dJC.updateDesignVariables(des);
std::string str;
if (startDes >= 0){
	char tmp[200];
	snprintf(tmp, 200,BASE_FOLDER "/../../results/n5/modL.txt");
	str = tmp;
    gsInfo << "Loading from " << str << "\n";
	des = loadVec(des.size(),str);
	dJC.updateDesignVariables(des);
}



  // gsInfo << dJC.getDvectors();
  // real_t eps = 0.5;
  //
	// std::srand((unsigned int) std::time(0));
  // gsVector<> ran;
	// ran.setRandom(des.size());
  // ran *= eps;

  gsInfo << "The domain is a "<< mp <<"\n";

  gsInfo << "evalObj() = " << evalObj2(des,quA,quB,mp,dJC);
  gsVector<> grad12 = gradObjImpl12(des,quA,quB,mp,dJC);
  gsVector<> grad22 = gradObjImpl22(des,quA,quB,mp,dJC);

  // gsMatrix<> disp(grad12.rows(),2);
  // disp << grad12, grad22;
  // gsInfo << "\n" << disp << "\n" ;

  // convergenceTestOfJacobian2(quA,quB,mp,dJC,impl);

  gsVector<> vec = evalVec(des,quA,quB,mp,dJC);

  convergenceTestOfJacobian3(quA,quB,mp,dJC,impl);


  exit(0);
  gsVector<> grad1 = gradObjImpl1(des,quA,quB,mp,dJC);
  gsVector<> grad2 = gradObjImpl2(des,quA,quB,mp,dJC);
  gsVector<> grad3 = gradObjImpl3(des,quA,quB,mp,dJC);

  gsInfo << "Difference in implementation 1 and 2 : " << (grad1-grad2).norm() << "\n";

  // gsMatrix<> disp(grad1.rows(),3);
  // disp << grad1, grad2, grad3;
  // gsInfo << "\n" << disp << "\n" ;

  convergenceTestOfJacobian(quA,quB,mp,dJC,impl);

  if (!output.empty()){
    gsInfo << "Output to " << output << "\n";
  }

  gsInfo << "\n\n==== DONE ====\n";
  return 0;
}
