#include <gismo.h>
#include "stateEquationAntenna.h"
using namespace gismo;


gsMatrix<> stateEquationAntenna::getDerivativeOfRhsZeroBC(index_t realOrImag){
	// gsInfo << "getDerivativeOfRhsZeroBC\n" << std::flush;
	// Method to generate system matrices and rhs
	// Input 0 for real part of matrix, and 1 for imaginary part of matrix

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

	geometryMap G = A.getMap(*mp);

	gsFunctionExpr<> zero("0.0",2);
	variable ffzero = A.getCoeff(zero,G);

	A.setIntegrationElements(dbasis);

  space u = A.getSpace(dbasis);
  u.setInterfaceCont(0);

	gsMultiBasis<> geom_basis(*mp);
  space v = A.getTestSpace(u,geom_basis);
	// Setup terms
	// Rhs
	variable Hiz_re = A.getCoeff(Hiz_real,G);
	variable Hiz_im = A.getCoeff(Hiz_imag,G);
	variable dHizdn_re = A.getCoeff(dHizdn_real,G);
	variable dHizdn_im = A.getCoeff(dHizdn_imag,G);

	variable dHizdx_re = A.getCoeff(dHizdx_real,G);
	variable dHizdx_im = A.getCoeff(dHizdx_imag,G);
	variable d2Hizdndx_re = A.getCoeff(d2Hizdndx_real,G);
	variable d2Hizdndx_im = A.getCoeff(d2Hizdndx_imag,G);
	// variable dHizdx_re = A.getCoeff();

	auto bnd_f_real = pde_eps_cr_inv*dHizdn_re + pde_bnd_const.real()*Hiz_re - pde_bnd_const.imag()*Hiz_im;
	auto bnd_f_imag = pde_eps_cr_inv*dHizdn_im + pde_bnd_const.imag()*Hiz_re + pde_bnd_const.real()*Hiz_im;

	auto bnd_dfdx_real = pde_eps_cr_inv*d2Hizdndx_re + pde_bnd_const.real()*dHizdx_re - pde_bnd_const.imag()*dHizdx_im;
	auto bnd_dfdx_imag = pde_eps_cr_inv*d2Hizdndx_im + pde_bnd_const.imag()*dHizdx_re + pde_bnd_const.real()*dHizdx_im;

	auto rhs_term_real = u*bnd_f_real*nv(G).norm();
	auto rhs_term_imag = u*bnd_f_imag*nv(G).norm();

	A.initSystem();

  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = A.getCoeff(x);
  variable fy = A.getCoeff(y);
  auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
  auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
  auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
  auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

  auto vxi = grad(v)*fjac(fx);
  auto veta = grad(v)*fjac(fy);

  auto detJ = j00*j11 - j10*j01;
	auto detJinv = 1/detJ.val();

	auto signOfDetJ = detJinv*meas(G);

	// FIXIT: veta here is due to the boundary we integrate over are east or west, TODO try to generalize...
	auto term_1x = veta*j01/nv(G).norm()*u.tr();
	auto term_2x = v*u.tr()*nv(G).norm();

	// gsInfo << "x terms\n" << std::flush;
	// A.assemble(term_1x);
	if (realOrImag == 0){
		// FIXIT: why including rhs here? TODO avoid this...
		A.assembleLhsRhsBc(bnd_f_real.val()*term_1x + bnd_dfdx_real.val()*term_2x,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
	} else {
		// FIXIT: why including rhs here? TODO avoid this...
		A.assembleLhsRhsBc(bnd_f_imag.val()*term_1x + bnd_dfdx_imag.val()*term_2x,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
	}

	gsMatrix<> xJac = A.matrix();

	auto term_1y = v*j11/nv(G).norm()*u.tr();

	// gsInfo << "y terms\n" << std::flush;
	A.initSystem();

	if (realOrImag == 0){
		// FIXIT: why including rhs here? TODO avoid this...
		A.assembleLhsRhsBc(bnd_f_real.val()*term_1y,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
	} else {
		// FIXIT: why including rhs here? TODO avoid this...
		A.assembleLhsRhsBc(bnd_f_imag.val()*term_1y,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
	}

	gsMatrix<> yJac = A.matrix();

	// gsInfo << "collect x and y terms\n";
	gsMatrix<> out;
	out.setZero(xJac.rows()*2,xJac.cols());

	// printMatSize(xJac,"xJac");
	// printMatSize(yJac,"yJac");
	// printMatSize(out,"out");

	out << xJac,
	 			 yJac;

  // delete f;
	// delete ms;
	// delete df_dx;
	// delete df_dy;
	return out;

}

gsMatrix<> stateEquationAntenna::getRhsZeroBC(index_t realOrImag){
	gsSparseMatrix<> mat;
	gsVector<> rhs;
	getTerm(realOrImag,mat,rhs);
	return rhs;
};
//
gsMatrix<> stateEquationAntenna::getDerivativeOfAu(index_t realOrImag, gsMultiPatch<> sol){
	// gsInfo << "f = " << *f << "\n" << std::flush;
	// gsInfo << "ms = " << *ms << "\n" << std::flush;

	gsFunctionExpr<> gN("0.0",2);
	gsFunctionExpr<> zero("0.0",2);
  //! [Boundary conditions]

	gsExprAssembler<> A(1,1);
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

	geometryMap G = A.getMap(*mp);

  gsExprEvaluator<> ev(A);
	A.setIntegrationElements(dbasis);

  space u = A.getSpace(dbasis);
  u.setInterfaceCont(0);
  u.addBc(bcInfoZero.get("Dirichlet"));

	gsMultiBasis<> geom_basis(*mp);

  space v = A.getTestSpace(u,geom_basis);
  // v.setInterfaceCont(0);

	// variable ff = ev.getVariable(f,G);
	// variable dff_dx = A.getCoeff(*df_dx, G);
	// variable dff_dy = A.getCoeff(*df_dy, G);

	variable solVar = A.getCoeff(sol);
	variable ffzero = A.getCoeff(zero,G);

	A.initSystem();

  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = A.getCoeff(x);
  variable fy = A.getCoeff(y);
  auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
  auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
  auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
  auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

  auto vxi = grad(v)*fjac(fx);
  auto veta = grad(v)*fjac(fy);

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

	variable geom = A.getCoeff(*mp);

	// auto G_x = fjac(fx).tr()*geom;
	// auto G_y = fjac(fy).tr()*geom;

	// gsInfo << "Initialize Solxi and Soleta\n";
  auto solxi = (fjac(solVar).tr()*jac(G).inv())*fjac(fx);
  auto soleta = (fjac(solVar).tr()*jac(G).inv())*fjac(fy);

	gsFunctionExpr<> ff3("y","-x",2);
	variable f3 = A.getCoeff(ff3);
	auto mat3 = fjac(f3);

  auto detJ = j00*j11 - j10*j01;
	auto detJinv = 1/detJ.val();

	auto signOfDetJ = detJinv*meas(G);

  auto d_detJ_dcx = signOfDetJ*(vxi*j11 - veta*j10) ;
  auto d_detJ_dcy = signOfDetJ*(veta*j00 - vxi*j01) ;

	auto igradsoltr = jac(G).inv().tr()*fjac(solVar);
	auto igradsol = fjac(solVar).tr()*jac(G).inv();

	auto igradsol0 = (igradsol*fjac(fx)).val();
	auto igradsol1 = (igradsol*fjac(fy)).val();

	variable fun_real = ev.getVariable(pde_eps_cr_fun_real);
	variable fun_imag = ev.getVariable(pde_eps_cr_fun_imag);

	auto laplace_term_real =	fun_real.val()*igrad(u,G)*igrad(u,G).tr()*meas(G);
	auto laplace_term_imag =	fun_imag.val()*igrad(u,G)*igrad(u,G).tr()*meas(G);

	// Helmholtz Term
	variable fun2_real = ev.getVariable(pde_mu_r_fun_real);
	real_t k0sq = pde_k0*pde_k0;

	auto helmholtz_term_real = -k0sq*fun2_real.val()*u*u.tr()*meas(G);

	// Boundary part
	auto bnd_term_real = pde_bnd_const.real()*u*u.tr()*nv(G).norm();
	auto bnd_term_imag = pde_bnd_const.imag()*u*u.tr()*nv(G).norm();

	// Rhs
	variable Hiz_re = A.getCoeff(Hiz_real,G);
	variable Hiz_im = A.getCoeff(Hiz_imag,G);
	variable dHizdn_re = A.getCoeff(dHizdn_real,G);
	variable dHizdn_im = A.getCoeff(dHizdn_imag,G);

	auto rhs_term_real = u*(pde_eps_cr_inv*dHizdn_re + pde_bnd_const.real()*Hiz_re - pde_bnd_const.imag()*Hiz_im)*nv(G).norm();
	auto rhs_term_imag = u*(pde_eps_cr_inv*dHizdn_im + pde_bnd_const.imag()*Hiz_re + pde_bnd_const.real()*Hiz_im)*nv(G).norm();


	// Laplace part
	//FIXIT check the signs on term 1 and 4
	auto term_1x = -d_detJ_dcx*igradsol*igrad(u,G).tr();
	auto term_2x = signOfDetJ*grad(v)*mat3*fjac(solVar)*(fjac(fy).tr()*igrad(u,G).tr());//(grad(v)*mat3)*soleta.val()*igrad(u,G).tr();
	auto term_3x = signOfDetJ*grad(v)*mat3*grad(u).tr()*igradsol1;//(grad(v)*mat3)*(jac(G).inv().tr()*fjac(solVar))*ueta.tr();;

	// Helmholtz part
	auto term_4x = -k0sq*d_detJ_dcx*solVar.val()*u.tr();

	// Bnd part
	auto term_5x = veta*j01/nv(G).norm()*solVar.val()*u.tr();

	// gsInfo << "x terms\n";

	A.initSystem();
	if (realOrImag == 0){
		A.assemble(fun_real.val()*(term_1x + term_2x + term_3x) + fun2_real.val()*term_4x);
		A.assembleLhsRhsBc(pde_bnd_const.real()*term_5x,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
	} else {
		A.assemble(fun_imag.val()*(term_1x + term_2x + term_3x));

		// A.assemble(laplace_term_imag) ;
		A.assembleLhsRhsBc(pde_bnd_const.imag()*term_5x, v*ffzero*nv(G).norm(), bcInfo.neumannSides());
	}
	gsMatrix<> xJac = A.matrix();

	A.initSystem();
	// A.initMatrix();

	// Laplace part
	auto term_1y = -d_detJ_dcy*igradsol*igrad(u,G).tr();
	auto term_2y = -signOfDetJ*grad(v)*mat3*fjac(solVar)*(fjac(fx).tr()*igrad(u,G).tr());//(grad(v)*mat3)*soleta.val()*igrad(u,G).tr();
	auto term_3y = -signOfDetJ*grad(v)*mat3*grad(u).tr()*igradsol0;//(grad(v)*mat3)*(jac(G).inv().tr()*fjac(solVar))*ueta.tr();;

	// Helmholtz part
	auto term_4y = -k0sq*d_detJ_dcy*solVar.val()*u.tr();

	// Bnd part
	auto term_5y = veta*j11/nv(G).norm()*solVar.val()*u.tr();

	// gsInfo << "y terms\n";
	if (realOrImag == 0){
		A.assemble(fun_real.val()*(term_1y + term_2y + term_3y) + fun2_real.val()*term_4y);
		A.assembleLhsRhsBc(pde_bnd_const.real()*term_5y,v*ffzero*nv(G).norm(), bcInfo.neumannSides());
	} else {
		A.assemble(fun_imag.val()*(term_1y + term_2y + term_3y));
		A.assembleLhsRhsBc(pde_bnd_const.imag()*term_5y,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
	}

	// gsInfo << "assembly done\n" << std::flush;

	gsMatrix<> yJac = A.matrix();

	// gsInfo << "collect x and y terms\n";
	gsMatrix<> out;
	out.setZero(xJac.rows()*2,xJac.cols());

	// printMatSize(xJac,"xJac");
	// printMatSize(yJac,"yJac");
	// printMatSize(out,"out");

	out << xJac,
	 			 yJac;

  // delete f;
	// delete ms;
	// delete df_dx;
	// delete df_dy;

	// gsInfo << "Return dKu \n";
	return out;

}
//
// gsMatrix<> stateEquationAntenna::getKu(gsMultiPatch<> sol){
// 	// gsInfo << "f = " << *f << "\n" << std::flush;
// 	// gsInfo << "ms = " << *ms << "\n" << std::flush;
//
// 	gsFunctionExpr<> gN("0.0",2);
// 	gsFunctionExpr<> zero("0.0",2);
//   //! [Boundary conditions]
//
// 	gsExprAssembler<> A(1,1);
//   geometryMap G = A.getMap(*mp);
//
// 	A.setIntegrationElements(dbasis);
//
//
//   space u = A.getSpace(dbasis);
//   u.setInterfaceCont(0);
//   u.addBc(bcInfoZero.get("Dirichlet"));
//
// 	variable solVar = A.getCoeff(sol);
//
// 	A.initSystem();
//
// 	A.assemble(igrad(u,G)*(jac(G).inv().tr()*fjac(solVar))*meas(G));
//
//   // delete f;
// 	// delete ms;
// 	// delete df_dx;
// 	// delete df_dy;
// 	gsInfo << "Return Ku \n";
// 	return A.rhs();
//
// }
//
gsMatrix<> stateEquationAntenna::getDerivativeOfU(){
	// gsInfo << "Get derivative of u\n";
	// gsInfo << "Get sol\n";
	gsMultiPatch<> u_real, u_imag;

	solve(u_real,u_imag);

	// gsInfo << "Get deriv of Ku\n";
	gsMatrix<> dK_realu_real = getDerivativeOfAu(0,u_real);
	gsMatrix<> dK_imagu_real = getDerivativeOfAu(1,u_real);
	gsMatrix<> dK_realu_imag = getDerivativeOfAu(0,u_imag);
	gsMatrix<> dK_imagu_imag = getDerivativeOfAu(1,u_imag);
	// gsInfo << "Get deriv of F\n";
	gsMatrix<> dF_real = getDerivativeOfRhsZeroBC(0);
	gsMatrix<> dF_imag = getDerivativeOfRhsZeroBC(1);

	gsMatrix<> dF(dF_real.rows(),dF_real.cols()*2);
	dF << dF_real, -dF_imag;

	gsMatrix<> dAu(dF_real.rows(),dF_real.cols()*2);
	dAu << dK_realu_real - dK_imagu_imag, -dK_imagu_real - dK_realu_imag;

	// gsInfo << "solve\n";
	gsMatrix<> du = solver.solve(dF.transpose() - dAu.transpose());

	// gsInfo << "return du\n";
	return du.transpose();
}

// Get the rhs of the sytem you need to solve to obtain dudc, without solving
gsMatrix<> stateEquationAntenna::getDerivativeWithoutSolving(){
	// gsInfo << "Get derivative of u\n";
	// gsInfo << "Get sol\n";
	gsMultiPatch<> u_real, u_imag;
	solve(u_real,u_imag);

	// gsInfo << "Get deriv of Ku\n";
	gsMatrix<> dK_realu_real = getDerivativeOfAu(0,u_real);
	gsMatrix<> dK_imagu_real = getDerivativeOfAu(1,u_real);
	gsMatrix<> dK_realu_imag = getDerivativeOfAu(0,u_imag);
	gsMatrix<> dK_imagu_imag = getDerivativeOfAu(1,u_imag);
	// gsInfo << "Get deriv of F\n";
	gsMatrix<> dF_real = getDerivativeOfRhsZeroBC(0);
	gsMatrix<> dF_imag = getDerivativeOfRhsZeroBC(1);

	gsMatrix<> dF(dF_real.rows(),dF_real.cols()*2);
	dF << dF_real, -dF_imag;

	gsMatrix<> dAu(dF_real.rows(),dF_real.cols()*2);
	dAu << dK_realu_real - dK_imagu_imag, -dK_imagu_real - dK_realu_imag;

	// Return matrix

	return dF.transpose() - dAu.transpose();

}

// Assumes the system is already factorized, e.g. if you have called getDerivativeWithoutSolving firts;
gsVector<> stateEquationAntenna::solveAdjoint(gsVector<> &rhs){
	return solver.solve(rhs);
}
//
// void stateEquationAntenna::plotMesh(gsMatrix<> solVector){
//
//
//   // gsMultiBasis<> dbasis(*mp);
//   // dbasis.setDegree(degree);
// 	//
//   // gsExprAssembler<> A(1,1);
//   // typedef gsExprAssembler<>::geometryMap geometryMap;
//   // typedef gsExprAssembler<>::variable    variable;
//   // typedef gsExprAssembler<>::space       space;
//   // typedef gsExprAssembler<>::solution    solution;
// 	//
//   // A.setIntegrationElements(dbasis);
//   // gsExprEvaluator<> ev(A);
// 	//
//   // geometryMap G = A.getMap(*mp);
// 	//
//   // space u = A.getSpace(dbasis);
//   // u.setInterfaceCont(0);
// 	//
//   // A.initSystem();
//   // solution u_sol = A.getSolution(u,solVector);
//
//   // mesh holds the control net of a geometry
//   // mesh is a set of vertices and lines (connections between vertices)
// 	for(index_t i = 0; i < mp->nBoxes(); i++){
//   	gsMesh<> mesh;
//   	mp->patch(i).controlNet(mesh);
// 		auto name = "mesh" + std::to_string( i );
//   	gsInfo << "Writing the control net to a paraview file: " <<  "\n" << "\n";
//   	gsWriteParaview(mesh, name);
// 	}
//
//   // gsInfo<<"Plotting in Paraview...\n";
//   // ev.options().setSwitch("plot.elements", false);
//   // ev.writeParaview( u_sol   , G, "solutionState");
//
// }
//
void stateEquationAntenna::printMatSize(gsMatrix<> mat, std::string name){
	gsInfo << "Size of " << name << ":\t (" << mat.rows() << ", " << mat.cols() << ")\n";
}
//
// void stateEquationAntenna::getFandMS(gsFunctionExpr<> *&fin, gsFunctionExpr<> *&msin, gsFunctionExpr<> *&df_dxin, gsFunctionExpr<> *&df_dyin){
// 	fin = &f;
// 	msin = &ms;
// 	df_dxin = &df_dx;
// 	df_dyin= &df_dy;
// }
//
// void stateEquationAntenna::getMSDerivatives(gsFunctionExpr<> *&dms_dxin, gsFunctionExpr<> *&dms_dyin){
// 	dms_dxin = &dms_dx;
// 	dms_dyin = &dms_dy;
// }
//
void stateEquationAntenna::plotSolution(gsMultiPatch<> &sol, std::string name){
	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);
  A.setIntegrationElements(dbasis);

  geometryMap G = A.getMap(*mp);

	variable u_sol = A.getCoeff(sol);

	gsInfo<<"Plotting " << name << " in Paraview...\n";
	// ev.options().setSwitch("plot.elements", true);
	ev.writeParaview( u_sol   , G, name);
}

void stateEquationAntenna::plotMagnitude(std::string name){
	gsMultiPatch<> u_real, u_imag;
	solve(u_real,u_imag);

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);
  A.setIntegrationElements(dbasis);

  geometryMap G = A.getMap(*mp);

	variable u_re = A.getCoeff(u_real);
	variable u_im = A.getCoeff(u_imag);

	gsInfo<<"Plotting " << name << " in Paraview...\n";
	// ev.options().setSwitch("plot.elements", true);
	ev.writeParaview( u_re*u_re + u_im*u_im   , G, name);
}

gsMultiPatch<> stateEquationAntenna::getPieceWiseFunctionOnPatch(index_t nBoxes, index_t patch, real_t val_on_patch, real_t val_elsewhere){
	gsKnotVector<> u_knots(0,1,0,1);
	gsKnotVector<> v_knots(0,1,0,1);


	gsTensorBSplineBasis<2> T_basis(u_knots,v_knots);

	gsMatrix<> cf_patch(1,1);
	cf_patch << val_on_patch;

	gsMatrix<> cf_elsewhere(1,1);
	cf_elsewhere << val_elsewhere;

	gsMultiPatch<> pw;

	for(index_t i = 0; i < nBoxes; i ++){
		if (i == patch){
  		memory::unique_ptr<gsGeometry<> > b = T_basis.makeGeometry(cf_patch);
			pw.addPatch(*b);
		} else {
  		memory::unique_ptr<gsGeometry<> > b = T_basis.makeGeometry(cf_elsewhere);
			pw.addPatch(*b);
		}
	}

	return pw;
}

void stateEquationAntenna::getTerm(index_t realOrImag, gsSparseMatrix<> &mat, gsVector<> &rhs){
	// Method to generate system matrices and rhs
	// Input 0 for real part of matrix, and 1 for imaginary part of matrix

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

	geometryMap G = A.getMap(*mp);

	A.setIntegrationElements(dbasis);

  space u = A.getSpace(dbasis);
  u.setInterfaceCont(0);

	// Setup terms

	// Laplace Term

	variable fun_real = ev.getVariable(pde_eps_cr_fun_real);
	variable fun_imag = ev.getVariable(pde_eps_cr_fun_imag);

	auto laplace_term_real =	fun_real.val()*igrad(u,G)*igrad(u,G).tr()*meas(G);
	auto laplace_term_imag =	fun_imag.val()*igrad(u,G)*igrad(u,G).tr()*meas(G);

	// Helmholtz Term

	variable fun2_real = ev.getVariable(pde_mu_r_fun_real);
	real_t k0sq = pde_k0*pde_k0;

	auto helmholtz_term_real = -k0sq*fun2_real.val()*u*u.tr()*meas(G);

	// Boundary part
	auto bnd_term_real = pde_bnd_const.real()*u*u.tr()*nv(G).norm();
	auto bnd_term_imag = pde_bnd_const.imag()*u*u.tr()*nv(G).norm();

	// Rhs
	variable Hiz_re = A.getCoeff(Hiz_real,G);
	variable Hiz_im = A.getCoeff(Hiz_imag,G);
	variable dHizdn_re = A.getCoeff(dHizdn_real,G);
	variable dHizdn_im = A.getCoeff(dHizdn_imag,G);

	auto rhs_term_real = u*(pde_eps_cr_inv*dHizdn_re + pde_bnd_const.real()*Hiz_re - pde_bnd_const.imag()*Hiz_im)*nv(G).norm();
	auto rhs_term_imag = u*(pde_eps_cr_inv*dHizdn_im + pde_bnd_const.imag()*Hiz_re + pde_bnd_const.real()*Hiz_im)*nv(G).norm();
	// auto rhs_term_real = 1.0*u*nv(G).norm();
	// auto rhs_term_imag = 1.0*u*nv(G).norm();

	A.initSystem();
	if (realOrImag == 0){
		// A.assemble(helmholtz_term_real);
		A.assemble(laplace_term_real + helmholtz_term_real);
		A.assembleLhsRhsBc(bnd_term_real,rhs_term_real, bcInfo.neumannSides());
		// A.assembleRhsBc(rhs_term_real, bcInfo.neumannSides());
	} else {
		A.assemble(laplace_term_imag) ;
		A.assembleLhsRhsBc(bnd_term_imag,rhs_term_imag, bcInfo.neumannSides());
		// A.assembleRhsBc(rhs_term_imag, bcInfo.neumannSides());
	}

	mat = A.matrix();
	rhs = A.rhs();
}

void stateEquationAntenna::getSystem(gsSparseMatrix<> &mat, gsVector<> &rhs){
	gsSparseMatrix<> matReal,matImag;
	gsVector<> rhsReal, rhsImag;

	// Get real parts
	getTerm(0,matReal,rhsReal);
	// gsInfo << "Save Ar \n" << std::flush;
	// std::ofstream file1("Ar.txt");
	// file1 << matReal.toDense();
	// file1.close();

	// Get imaginary parts
	getTerm(1,matImag,rhsImag);
	// gsInfo << "Save Ai \n" << std::flush;
	// std::ofstream file("Ai.txt");
	// file << matImag.toDense();
	// file.close();

	index_t rows = matReal.rows();
	index_t cols = matReal.cols();
	// gsInfo << "DoFs :" << rows*2 << "\n";
	gsSparseMatrix<> A(rows*2,cols*2);
	A.reserve(matReal.nonZeros()*2 + matImag.nonZeros()*2);
	for(index_t c = 0; c < cols; ++c)
	{
	    A.startVec(c); // Important: Must be called once for each column before inserting!

	    for(gsSparseMatrix<>::InnerIterator itR(matReal, c); itR; ++itR){
	         A.insertBack(itR.row(), c) = itR.value();
				 }
	    for(gsSparseMatrix<>::InnerIterator itI(matImag, c); itI; ++itI){
	         A.insertBack(itI.row()+rows, c) = -itI.value();
				 }
	}
	// gsInfo << " next loop...sss\n" << std::flush;
	for(index_t c = 0; c < cols; ++c)
	{
	    A.startVec(c+cols); // Important: Must be called once for each column before inserting!

	    for(gsSparseMatrix<>::InnerIterator itI(matImag, c); itI; ++itI){
	         A.insertBack(itI.row(), c + cols) = -itI.value();
				 }
	    for(gsSparseMatrix<>::InnerIterator itR(matReal, c); itR; ++itR){
	         A.insertBack(itR.row() + rows, c + cols) = -itR.value();
				 }
	}
	// gsInfo << "Finalize A \n" << std::flush;
	A.finalize();

	mat = A;
	rhs.setZero(cols*2);
	rhs << rhsReal, -rhsImag;

}

void stateEquationAntenna::solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		// gsInfo << "Save A \n" << std::flush;
		// std::ofstream file("A.txt");
		// file << mat.toDense();
		// file.close();

		gsInfo << "." << std::flush;
		solver.compute(mat);
		solVector = solver.solve(rhs);


		// Plot
		gsExprAssembler<> A(1,1);
		gsExprEvaluator<> ev(A);

		geometryMap G = A.getMap(*mp);

		A.setIntegrationElements(dbasis);

		space u = A.getSpace(dbasis);
		u.setInterfaceCont(0);
		// u.addBc(bcInfo.get("Dirichlet"));

		A.initSystem();

		gsMatrix<> solVector_Real = solVector.block(0,0,mat.cols()/2,1);
		gsMatrix<> solVector_Imag = solVector.block(mat.cols()/2,0,mat.cols()/2,1);

		solution u_sol_real = A.getSolution(u,solVector_Real);

		u_sol_real.extract(u_real);

		solution u_sol_imag = A.getSolution(u,solVector_Imag);

		u_sol_imag.extract(u_imag);
}

gsMatrix<> stateEquationAntenna::getU(index_t realOrImag){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		gsInfo << "." << std::flush;
		solver.compute(mat);
		solVector = solver.solve(rhs);

		index_t n = mat.cols()/2;

		gsMatrix<> U_Real = solVector.block(0,0,n,1);
		gsMatrix<> U_Imag = solVector.block(n,0,n,1);

		if (realOrImag == 0){
			return U_Real;
		} else {
			return U_Imag;
		}
}

gsMatrix<> stateEquationAntenna::getU(){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		gsInfo << "." << std::flush;
		solver.compute(mat);
		solVector = solver.solve(rhs);

		return solVector;
}

gsMatrix<> stateEquationAntenna::getAu(index_t realOrImag, gsMatrix<> U){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);

		index_t n = mat.cols()/2;

		if (realOrImag == 0){
			return mat.block(0,0,n,n)*U;
		} else {
			return mat.block(n,0,n,n)*U;
		}
}

void stateEquationAntenna::assembleAndSolve(){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);
		solver.compute(mat);
		gsInfo << "SOLVE:\n" << std::flush;
		solVector = solver.solve(rhs);

		// Plot
		gsExprAssembler<> A(1,1);
		gsExprEvaluator<> ev(A);

		geometryMap G = A.getMap(*mp);

		A.setIntegrationElements(dbasis);

		space u = A.getSpace(dbasis);
		u.setInterfaceCont(0);

		gsMatrix<> solVector_Real = solVector.block(0,0,mat.cols()/2,1);
		gsMatrix<> solVector_Imag = solVector.block(mat.cols()/2,0,mat.cols()/2,1);

		gsInfo << solVector_Real;

    solution u_sol_real = A.getSolution(u, solVector_Real);

    ev.options().setSwitch("plot.elements", true);
		gsInfo << "Plotting Solution \n" << std::flush;
    ev.writeParaview( u_sol_real   , G, "solutionHelmholtz");
		gsInfo << "Done Plotting \n" << std::flush;

	}

void stateEquationAntenna::printConstants(){
		gsInfo << "\n\n CONSTANTS: \n\n";

	  gsInfo << "r_t = " << pde_r_t << "\n";
    gsInfo << "k0 = " << pde_k0 << "\n";
    gsInfo << "f = " << pde_f << "\n";
    gsInfo << "f_r = " << pde_f_r << "\n";
	  gsInfo << "L_f = " << pde_L_f << "\n";
		gsInfo << "sigma = " << pde_sigma << "\n";
		gsInfo << "omega = " << pde_omega << "\n";
		gsInfo << "eps_rs_real = " << pde_eps_rs_real << "\n";
		gsInfo << "eps_rs_imag = " << pde_eps_rs_imag << "\n";
		gsInfo << "eps_r = " << pde_eps_r << "\n";
    gsInfo << "mu0 = " << pde_mu0 << "\n";
    gsInfo << "c = " << pde_c << "\n";
    gsInfo << "eps0 = " << pde_eps0 << "\n";
		gsInfo << "eps_crs_real = " << pde_eps_crs_real << "\n";
		gsInfo << "eps_crs_imag = " << pde_eps_crs_imag << "\n";
		gsInfo << "eps_cr = " << pde_eps_cr << "\n";

		gsInfo << "mu_r = " << pde_mu_r << "\n";
		gsInfo << "mu_rs = " << pde_mu_rs << "\n";
		gsInfo << "mu_cr = " << pde_mu_cr << "\n";

    gsInfo << "term = " << sqrt(pde_eps0*pde_mu0) << "\n";
    gsInfo << "expon = " << pde_k0*sqrt(pde_eps_cr*pde_mu_r) << "\n";
    gsInfo << "pde_eps_crs_real = " << pde_eps_crs_real << "\n";
    gsInfo << "pde_eps_crs_imag = " << pde_eps_crs_imag << "\n";
		gsInfo << "\n\n ========== \n\n";

}
