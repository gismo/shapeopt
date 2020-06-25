#include <gismo.h>
#include "gsStateEquationAntenna.h"
using namespace gismo;


// Overloaded methods
gsMatrix<> gsStateEquationAntenna::getDerivativeOfRhsZeroBC(index_t realOrImag){
	// gsInfo << "getDerivativeOfRhsZeroBC\n" << std::flush;
	// Method to generate system matrices and rhs
	// Input 0 for real part of matrix, and 1 for imaginary part of matrix

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

    A.options().setReal("quA",m_quA);
    A.options().setInt("quB",m_quB);

	geometryMap G = A.getMap(*m_mp);

	gsFunctionExpr<> zero("0.0",2);
	variable ffzero = A.getCoeff(zero,G);

	A.setIntegrationElements(dbasis);

  space u = A.getSpace(dbasis);
  u.setInterfaceCont(0);

	gsMultiBasis<> geom_basis(*m_mp);
  space v = A.getTestSpace(u,geom_basis);
	// Setup terms
	// Rhs
	variable Hiz_re = A.getCoeff(*Hiz_real,G);
	variable Hiz_im = A.getCoeff(*Hiz_imag,G);
	variable dHizdn_re = A.getCoeff(*dHizdn_real,G);
	variable dHizdn_im = A.getCoeff(*dHizdn_imag,G);

	variable dHizdx_re = A.getCoeff(*dHizdx_real,G);
	variable dHizdx_im = A.getCoeff(*dHizdx_imag,G);
	variable d2Hizdndx_re = A.getCoeff(*d2Hizdndx_real,G);
	variable d2Hizdndx_im = A.getCoeff(*d2Hizdndx_imag,G);
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
	auto term_1x = vxi*j00/nv(G).norm()*u.tr();
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

	auto term_1y = v*j10/nv(G).norm()*u.tr();

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

void gsStateEquationAntenna::getTerm(index_t realOrImag, gsSparseMatrix<> &mat, gsVector<> &rhs){
	// Method to generate system matrices and rhs
	// Input 0 for real part of matrix, and 1 for imaginary part of matrix

	gsExprAssembler<> A(1,1);
	gsExprEvaluator<> ev(A);

    A.options().setReal("quA",m_quA);
    A.options().setInt("quB",m_quB);

    geometryMap G = A.getMap(*m_mp);

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
	variable Hiz_re = A.getCoeff(*Hiz_real,G);
	variable Hiz_im = A.getCoeff(*Hiz_imag,G);
	variable dHizdn_re = A.getCoeff(*dHizdn_real,G);
	variable dHizdn_im = A.getCoeff(*dHizdn_imag,G);

	auto rhs_term_real = u*(pde_eps_cr_inv*dHizdn_re + pde_bnd_const.real()*Hiz_re - pde_bnd_const.imag()*Hiz_im)*nv(G).norm();
	auto rhs_term_imag = u*(pde_eps_cr_inv*dHizdn_im + pde_bnd_const.imag()*Hiz_re + pde_bnd_const.real()*Hiz_im)*nv(G).norm();
	// auto rhs_term_real = 1.0*u*nv(G).norm();
	// auto rhs_term_imag = 1.0*u*nv(G).norm();

	A.initSystem();
	if (realOrImag == 0){
		// A.assemble(helmholtz_term_real);
		A.assemble(laplace_term_real + helmholtz_term_real);
		A.assembleLhsRhsBc(bnd_term_real,rhs_term_real, bcInfo.neumannSides());
		//gsInfo << "A.matrix().norm()" << A.matrix() << "\n";
		// A.assembleRhsBc(rhs_term_real, bcInfo.neumannSides());
	} else {
		A.assemble(laplace_term_imag) ;
		A.assembleLhsRhsBc(bnd_term_imag,rhs_term_imag, bcInfo.neumannSides());
		// A.assembleRhsBc(rhs_term_imag, bcInfo.neumannSides());
	}

	mat = A.matrix();
	rhs = A.rhs();
}

void gsStateEquationAntenna::printConstants(){
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

// Local methods

gsMatrix<> gsStateEquationAntenna::getDerivativeOfAuPart2(index_t realOrImag, gsMultiPatch<> sol){
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

    A.options().setReal("quA",m_quA);
    A.options().setInt("quB",m_quB);

    geometryMap G = A.getMap(*m_mp);

    gsExprEvaluator<> ev(A);
    A.setIntegrationElements(dbasis);

    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);
    u.addBc(bcInfoZero.get("Dirichlet"));

    gsMultiBasis<> geom_basis(*m_mp);

    space v = A.getTestSpace(u,geom_basis,2);
    // v.setInterfaceCont(0);

    // variable ff = ev.getVariable(f,G);
    // variable dff_dx = A.getCoeff(*df_dx, G);
    // variable dff_dy = A.getCoeff(*df_dy, G);

    variable solVar = A.getCoeff(sol);
    variable ffzero = A.getCoeff(zero,G);

    A.initSystem();

    auto igradsoltr = jac(G).inv().tr()*fjac(solVar);
    auto igradsol = fjac(solVar).tr()*jac(G).inv();

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
    variable Hiz_re = A.getCoeff(*Hiz_real,G);
    variable Hiz_im = A.getCoeff(*Hiz_imag,G);
    variable dHizdn_re = A.getCoeff(*dHizdn_real,G);
    variable dHizdn_im = A.getCoeff(*dHizdn_imag,G);

    auto rhs_term_real = u*(pde_eps_cr_inv*dHizdn_re + pde_bnd_const.real()*Hiz_re - pde_bnd_const.imag()*Hiz_im)*nv(G).norm();
    auto rhs_term_imag = u*(pde_eps_cr_inv*dHizdn_im + pde_bnd_const.imag()*Hiz_re + pde_bnd_const.real()*Hiz_im)*nv(G).norm();

    // We use -meas(G) to ensure the same sign as meas(G).. (so it also works for negative determinant)
    // For easier reading, replace with - jac(G).det().sgn()*jac(G).det()
    // auto term_1 = meas(G)*matrix_by_space(jac(G).inv(),jac(v)).trace()*igradsol*igrad(u,G).tr();
    auto term_1 = matrix_by_space(jac(G).inv(),jac(v)).trace()*fjac(solVar).tr()*jac(G).inv()*igrad(u,G).tr()*meas(G);

    // auto term_2y = -signOfDetJ*grad(v)*mat3*fjac(solVar)*(fjac(fx).tr()*igrad(u,G).tr());//(grad(v)*mat3)*soleta.val()*igrad(u,G).tr();
    // auto term_2 =  -collapse(matrix_by_space_tr(jac(G).inv(),jac(v)),(jac(G).inv().tr()*fjac(solVar)))*igrad(u,G).tr()*meas(G);
    // auto term_2 =  -collapse(fjac(solVar).tr(),matrix_by_space(jac(G).inv(),jac(v))*jac(G).inv())*igrad(u,G).tr()*meas(G);
    auto term_2 = - collapse(fjac(solVar).tr(),matrix_by_space(jac(G).inv(),jac(v))*jac(G).inv()) * igrad(u,G).tr()*meas(G);
    // auto term_3 =  -collapse(igradsol,matrix_by_space_tr(jac(G).inv(),jac(v)))*igrad(u,G).tr()*meas(G);
    auto term_3 = - collapse(fjac(solVar).tr()*jac(G).inv(),matrix_by_space_tr(jac(G).inv(),jac(v)))*igrad(u,G).tr()*meas(G);
    // A.assemble(term_3);
    // auto term_2 = (jac(v)*fjac(solVar)).tr()*igrad(u,G).tr();
    // auto term_3 = (matrix_by_space(jac(G).inv(),jac(v))).tr()*igradsol*igrad(u,G).tr()*meas(G);

    // auto term_4x = -k0sq*d_detJ_dcx*solVar.val()*u.tr();
    auto term_4 = -k0sq*meas(G)*matrix_by_space(jac(G).inv(),jac(v)).trace()*solVar.val()*u.tr();

    // A.assemble(term_2);
    // gsMatrix<> tmp = A.matrix();
    // writeToFile(tmp,"../results/A2.txt");
    // exit(0);

    A.initSystem();
    if (realOrImag == 0){
        A.assemble(fun_real.val()*(term_3 + term_2 + term_1) + fun2_real.val()*term_4);
        // A.assembleLhsRhsBc(pde_bnd_const.real()*term_5x,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
    } else {
        A.assemble(fun_imag.val()*(term_3 + term_2 + term_1));
        // A.assemble(laplace_term_imag) ;
        // A.assembleLhsRhsBc(pde_bnd_const.imag()*term_5x, v*ffzero*nv(G).norm(), bcInfo.neumannSides());
    }

    // gsInfo << "\nsize (" << A.matrix().rows() << ", " << A.matrix().cols() << ")\n";

    // gsInfo << "collect x and y terms\n";
    // printMatSize(xJac,"xJac");
    // printMatSize(yJac,"yJac");
    // printMatSize(out,"out");

    // out << xJac,
    // yJac;

    // delete f;
    // delete ms;
    // delete df_dx;
    // delete df_dy;

    // gsInfo << "Return dKu \n";
    return A.matrix();

}

gsMatrix<> gsStateEquationAntenna::getDerivativeOfAuPart1(index_t realOrImag, gsMultiPatch<> sol){
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

    geometryMap G = A.getMap(*m_mp);

    gsExprEvaluator<> ev(A);
    A.setIntegrationElements(dbasis);

    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);
    u.addBc(bcInfoZero.get("Dirichlet"));

    gsMultiBasis<> geom_basis(*m_mp);

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

    variable geom = A.getCoeff(*m_mp);

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
    variable Hiz_re = A.getCoeff(*Hiz_real,G);
    variable Hiz_im = A.getCoeff(*Hiz_imag,G);
    variable dHizdn_re = A.getCoeff(*dHizdn_real,G);
    variable dHizdn_im = A.getCoeff(*dHizdn_imag,G);

    auto rhs_term_real = u*(pde_eps_cr_inv*dHizdn_re + pde_bnd_const.real()*Hiz_re - pde_bnd_const.imag()*Hiz_im)*nv(G).norm();
    auto rhs_term_imag = u*(pde_eps_cr_inv*dHizdn_im + pde_bnd_const.imag()*Hiz_re + pde_bnd_const.real()*Hiz_im)*nv(G).norm();


    // Laplace part
    //FIXIT check the signs on term 1 and 4
    // auto term_1x = -d_detJ_dcx*igradsol*igrad(u,G).tr();
    // auto term_2x = signOfDetJ*grad(v)*mat3*fjac(solVar)*(fjac(fy).tr()*igrad(u,G).tr());//(grad(v)*mat3)*soleta.val()*igrad(u,G).tr();
    // auto term_3x = signOfDetJ*grad(v)*mat3*grad(u).tr()*igradsol1;//(grad(v)*mat3)*(jac(G).inv().tr()*fjac(solVar))*ueta.tr();;

    // Helmholtz part
    // auto term_4x = -k0sq*d_detJ_dcx*solVar.val()*u.tr();

    // Bnd part
    // auto term_5x = veta*j01/nv(G).norm()*solVar.val()*u.tr();
	auto term_5x = vxi*j00/nv(G).norm()*solVar.val()*u.tr();
    // gsInfo << "x terms\n";

    A.initSystem();
    if (realOrImag == 0){
        // A.assemble(fun_real.val()*(term_2x));
        A.assembleLhsRhsBc(pde_bnd_const.real()*term_5x,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
    } else {
        // A.assemble(fun_imag.val()*(term_2x));

        // A.assemble(laplace_term_imag) ;
        A.assembleLhsRhsBc(pde_bnd_const.imag()*term_5x, v*ffzero*nv(G).norm(), bcInfo.neumannSides());
    }
    gsMatrix<> xJac = A.matrix();

    A.initSystem();
    // A.initMatrix();

    // Laplace part
    // auto term_1y = -d_detJ_dcy*igradsol*igrad(u,G).tr();
    // auto term_2y = -signOfDetJ*grad(v)*mat3*fjac(solVar)*(fjac(fx).tr()*igrad(u,G).tr());//(grad(v)*mat3)*soleta.val()*igrad(u,G).tr();
    // auto term_3y = -signOfDetJ*grad(v)*mat3*grad(u).tr()*igradsol0;//(grad(v)*mat3)*(jac(G).inv().tr()*fjac(solVar))*ueta.tr();;

    // Helmholtz part
    // auto term_4y = -k0sq*d_detJ_dcy*solVar.val()*u.tr();

    // Bnd part
    // auto term_5y = veta*j11/nv(G).norm()*solVar.val()*u.tr();
	auto term_5y = vxi*j10/nv(G).norm()*solVar.val()*u.tr();

    // gsInfo << "y terms\n";
    if (realOrImag == 0){
        // A.assemble(fun_real.val()*(term_2y));
        A.assembleLhsRhsBc(pde_bnd_const.real()*term_5y,v*ffzero*nv(G).norm(), bcInfo.neumannSides());
    } else {
        // A.assemble(fun_imag.val()*(term_2y));
        A.assembleLhsRhsBc(pde_bnd_const.imag()*term_5y,v*ffzero*nv(G).norm(),bcInfo.neumannSides());
    }

    // gsInfo << "assembly done\n" << std::flush;

    gsMatrix<> yJac = A.matrix();

    // A.initSystem();
    // A.assemble(term_2x);
    // xJac = A.matrix();
    // A.initSystem();
    // A.assemble(term_2y);
    // yJac = A.matrix();

    // gsInfo << "collect x and y terms\n";
    gsMatrix<> out;
    out.setZero(xJac.rows()*2,xJac.cols());
    out << xJac,yJac;
    // writeToFile(out,"../results/A1.txt");
    // exit(0);
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

void gsStateEquationAntenna::assembleAndSolve(){
		gsSparseMatrix<> mat;
		gsVector<> rhs;
		getSystem(mat,rhs);
		solver.compute(mat);
		gsInfo << "SOLVE:\n" << std::flush;
		solVector = solver.solve(rhs);

		// Plot
		gsExprAssembler<> A(1,1);
		gsExprEvaluator<> ev(A);

		geometryMap G = A.getMap(*m_mp);

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

// Old methods

/* gsMatrix<> gsStateEquationAntenna::getKu(gsMultiPatch<> sol){
 	// gsInfo << "f = " << *f << "\n" << std::flush;
 	// gsInfo << "ms = " << *ms << "\n" << std::flush;

 	gsFunctionExpr<> gN("0.0",2);
 	gsFunctionExpr<> zero("0.0",2);
   //! [Boundary conditions]

 	gsExprAssembler<> A(1,1);
   geometryMap G = A.getMap(*m_mp);

 	A.setIntegrationElements(dbasis);


   space u = A.getSpace(dbasis);
   u.setInterfaceCont(0);
   u.addBc(bcInfoZero.get("Dirichlet"));

 	variable solVar = A.getCoeff(sol);

 	A.initSystem();

 	A.assemble(igrad(u,G)*(jac(G).inv().tr()*fjac(solVar))*meas(G));

   // delete f;
 	// delete ms;
 	// delete df_dx;
 	// delete df_dy;
 	gsInfo << "Return Ku \n";
 	return A.rhs();

 }
*/

/* void gsStateEquationAntenna::plotMesh(gsMatrix<> solVector){


   // gsMultiBasis<> dbasis(*m_mp);
   // dbasis.setDegree(degree);
 	//
   // gsExprAssembler<> A(1,1);
   // typedef gsExprAssembler<>::geometryMap geometryMap;
   // typedef gsExprAssembler<>::variable    variable;
   // typedef gsExprAssembler<>::space       space;
   // typedef gsExprAssembler<>::solution    solution;
 	//
   // A.setIntegrationElements(dbasis);
   // gsExprEvaluator<> ev(A);
 	//
   // geometryMap G = A.getMap(*m_mp);
 	//
   // space u = A.getSpace(dbasis);
   // u.setInterfaceCont(0);
 	//
   // A.initSystem();
   // solution u_sol = A.getSolution(u,solVector);

   // mesh holds the control net of a geometry
   // mesh is a set of vertices and lines (connections between vertices)
 	for(index_t i = 0; i < m_mp->nBoxes(); i++){
   	gsMesh<> mesh;
   	m_mp->patch(i).controlNet(mesh);
 		auto name = "mesh" + std::to_string( i );
   	gsInfo << "Writing the control net to a paraview file: " <<  "\n" << "\n";
   	gsWriteParaview(mesh, name);
 	}

   // gsInfo<<"Plotting in Paraview...\n";
   // ev.options().setSwitch("plot.elements", false);
   // ev.writeParaview( u_sol   , G, "solutionState");

 }
*/

/* void gsStateEquationAntenna::getFandMS(gsFunctionExpr<> *&fin, gsFunctionExpr<> *&msin, gsFunctionExpr<> *&df_dxin, gsFunctionExpr<> *&df_dyin){
 	fin = &f;
 	msin = &ms;
 	df_dxin = &df_dx;
 	df_dyin= &df_dy;
 }
*/

/*
 void gsStateEquationAntenna::getMSDerivatives(gsFunctionExpr<> *&dms_dxin, gsFunctionExpr<> *&dms_dyin){
 	dms_dxin = &dms_dx;
 	dms_dyin = &dms_dy;
 }
*/


