#ifndef GSSTATEEQUATIONANTENNA_H
#define GSSTATEEQUATIONANTENNA_H

#include <math.h>
#include <complex>
#include <sstream>
#include <string>
#include "gsStateEquation.h"
using namespace gismo;

class gsStateEquationAntenna : public gsStateEquation{
public:
    gsStateEquationAntenna(memory::shared_ptr<gsMultiPatch<>> mpin, index_t numRefine): 
        gsStateEquation(mpin, numRefine), zero("0.0",2)
    {
        pde_f       = pde_f/pde_L_f;
        pde_omega   = 2*M_PI*pde_f;
        pde_eps_rs_real  = -20.198689873114;
        pde_eps_rs_imag = 1.38170639237;
        pde_eps_crs_real = pde_eps_rs_real;
        pde_eps_crs_imag = pde_eps_rs_imag - pde_sigma/(pde_omega*pde_L_f*pde_eps0);
        pde_k0      = 2*M_PI*pde_f*sqrt(pde_eps0*pde_mu0);

        pde_eps_rs = pde_eps_rs_real + 1i*pde_eps_rs_imag;
        pde_eps_crs = pde_eps_crs_real + 1i*pde_eps_crs_imag;

        pde_bnd_const = 1.0/pde_eps_cr*(1i*pde_k0 + 1.0/(2*pde_r_t));

        pde_eps_crs_inv = 1.0/pde_eps_crs;
        pde_eps_cr_inv = 1.0/pde_eps_cr;

        //FIXIT Antenna patch is hardcoded to 3.!!!
        pde_eps_cr_fun_real = getPieceWiseFunctionOnPatch(m_mp->nBoxes(),3,pde_eps_crs_inv.real(), pde_eps_cr_inv);
        pde_eps_cr_fun_imag = getPieceWiseFunctionOnPatch(m_mp->nBoxes(),3,pde_eps_crs_inv.imag(), 0.0);
        // pde_eps_cr_fun_real = getPieceWiseFunctionOnPatch(m_mp->nBoxes(),3,pde_eps_cr_inv, pde_eps_cr_inv);
        // pde_eps_cr_fun_imag = getPieceWiseFunctionOnPatch(m_mp->nBoxes(),3,0.0, 0.0);
        pde_mu_r_fun_real = getPieceWiseFunctionOnPatch(m_mp->nBoxes(),3,pde_mu_rs, pde_mu_r);
        // plotSolution(pde_eps_cr_fun_real,"pde_eps_cr_fun_real");
        // plotSolution(pde_eps_cr_fun_imag,"pde_eps_cr_fun_imag");

        // printConstants();

        // Change domain size
        // for (index_t i = 0; i < m_mp->nPatches(); i++){
        //    gsMatrix<> cc = m_mp->patch(i).coefs();
        //    m_mp->patch(i).setCoefs(pde_r_t/8.0*cc);
        // }

        char tmp[200];
        snprintf(tmp, 200, "(-%f*sqrt(%f*%f))", pde_k0, pde_eps_cr, pde_mu_r);
        std::string k_str = tmp;

        snprintf(tmp, 200, "cos(%s*x)", k_str.c_str());
        std::string Hiz_real_str = tmp;
        // gsInfo << "\n Hiz real: " << Hiz_real_str << "\n";
        Hiz_real = memory::make_shared(new gsFunctionExpr<>(Hiz_real_str,2));

        snprintf(tmp, 200, "sin(%s*x)", k_str.c_str());
        std::string Hiz_imag_str = tmp;
        // gsInfo << "\n Hiz imag: " << Hiz_imag_str << "\n";
        Hiz_imag = memory::make_shared(new gsFunctionExpr<>(Hiz_imag_str,2));

        snprintf(tmp, 200, "(-%s*sin(%s*x))", k_str.c_str(), k_str.c_str());
        std::string dHizdx_real_str = tmp;
        // gsInfo << "\n dHizdx real: " << dHizdx_real_str << "\n";
        dHizdx_real = memory::make_shared(new gsFunctionExpr<>(dHizdx_real_str,2));

        snprintf(tmp, 200, "(%s*cos(%s*x))", k_str.c_str(), k_str.c_str());
        std::string dHizdx_imag_str = tmp;
        // gsInfo << "\n dHizdx imag: " << dHizdx_imag_str << "\n";
        dHizdx_imag = memory::make_shared(new gsFunctionExpr<>(dHizdx_imag_str,2));

        snprintf(tmp, 200, "x/%f*%s", pde_r_t, dHizdx_real_str.c_str());
        std::string dHizdn_real_str = tmp;
        // gsInfo << "\n dHizdn real: " << dHizdn_real_str << "\n";
        dHizdn_real = memory::make_shared(new gsFunctionExpr<>(dHizdn_real_str,2));

        snprintf(tmp, 200, "x/%f*%s", pde_r_t, dHizdx_imag_str.c_str());
        std::string dHizdn_imag_str = tmp;
        // gsInfo << "\n dHizdn imag: " << dHizdn_imag_str << "\n";
        dHizdn_imag = memory::make_shared(new gsFunctionExpr<>(dHizdn_imag_str,2));

        snprintf(tmp, 200, "1/%f*%s - x/%f*%s^2*%s", pde_r_t, dHizdx_real_str.c_str(), pde_r_t, k_str.c_str(), Hiz_real_str.c_str());
        std::string d2Hizdndx_real_str = tmp;
        // gsInfo << "\n d2Hizdndx real: " << d2Hizdndx_real_str << "\n";
        d2Hizdndx_real = memory::make_shared(new gsFunctionExpr<>(d2Hizdndx_real_str,2));

        snprintf(tmp, 200, "1/%f*%s - x/%f*%s^2*%s", pde_r_t, dHizdx_imag_str.c_str(), pde_r_t, k_str.c_str(), Hiz_imag_str.c_str());
        std::string d2Hizdndx_imag_str = tmp;
        // gsInfo << "\n dHizdn imag: " << d2Hizdndx_imag_str << "\n";
        d2Hizdndx_imag = memory::make_shared(new gsFunctionExpr<>(d2Hizdndx_imag_str,2));

        bcInfo.addCondition(4, boundary::south, condition_type::neumann, &zero);
        // bcInfoZero.addCondition(4, boundary::east, condition_type::neumann, &zero); //zero bcs
        // bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &Hiz_real);
        // bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, &Hiz_real);
        // bcInfo.addCondition(2, boundary::west,  condition_type::dirichlet, &Hiz_real);
        // bcInfo.addCondition(4, boundary::east, condition_type::dirichlet, &Hiz_real);
    };

// Overloaded methods
    gsMatrix<> getDerivativeOfAu(index_t realOrImag, gsMultiPatch<> sol){
        return getDerivativeOfAuPart1(realOrImag,sol) + getDerivativeOfAuPart2(realOrImag,sol);
    };

    gsMatrix<> getDerivativeOfRhsZeroBC(index_t realOrImag);

    void getTerm(index_t realOrImag, gsSparseMatrix<> &mat, gsVector<> &rhs);

    void printConstants();

// Local methods
    
    // For implenting new gradients
    gsMatrix<> getDerivativeOfAuPart1(index_t realOrImag, gsMultiPatch<> sol); 

    gsMatrix<> getDerivativeOfAuPart2(index_t realOrImag, gsMultiPatch<> sol);

    void assembleAndSolve();

// Old methods

    // gsMatrix<> getKu(gsMultiPatch<> sol);
    // void plotMesh(gsMatrix<> solVector);
    // void getFandMS(gsFunctionExpr<> *&f, gsFunctionExpr<> *&ms, gsFunctionExpr<> *&df_dx, gsFunctionExpr<> *&df_dy);
    // void getMSDerivatives(gsFunctionExpr<> *&dms_dx, gsFunctionExpr<> *&dms_dy);

public:
    // Boundary conditions
    gsFunctionExpr<> zero;

    gsFunctionExpr<>::Ptr Hiz_real;
    gsFunctionExpr<>::Ptr Hiz_imag;

    gsFunctionExpr<>::Ptr dHizdx_real;
    gsFunctionExpr<>::Ptr dHizdx_imag;

    gsFunctionExpr<>::Ptr dHizdn_real;
    gsFunctionExpr<>::Ptr dHizdn_imag;
    gsFunctionExpr<>::Ptr d2Hizdndx_real;
    gsFunctionExpr<>::Ptr d2Hizdndx_imag;

    gsBoundaryConditions<> bcInfo;
    gsBoundaryConditions<> bcInfoZero;

    // Physical parameters
    real_t pde_r_t     = 4.0;                               // ok
    real_t pde_f       = 4e14;                              // ok
    real_t pde_f_r     = 1.15e8;                            // ok
    real_t pde_L_f     = pde_f/pde_f_r;                             // ok, but unsure where to use
    real_t pde_sigma   = 1e6;                               // ok
    real_t pde_omega   = 2*M_PI*pde_f;                          // ok
    real_t pde_eps_rs_real  = -20.198689873114;
    real_t pde_eps_rs_imag = 1.38170639237;                 // ok
    real_t pde_eps_r   = 1.0;                               // ok
    real_t pde_mu0     = 4*M_PI*1e-7;                       // From wikipedia, [H/m]
    real_t pde_c       = 299792458;                         // Speed of light, wikipedia, [m/s]
    real_t pde_eps0    = 1/(pde_mu0*pde_c*pde_c);                       // From wikipedia, [F/m]
    real_t pde_eps_cr  = pde_eps_r;                             // ok, found similar value on wikipedia -> relative permittivity
    real_t pde_eps_crs_real = pde_eps_rs_real;                  // ok, found similar value at refractiveindex.info
    real_t pde_eps_crs_imag = pde_eps_rs_imag;                  // ok, found similar value at refractiveindex.info
    real_t pde_mu_r    = 1.0;                               // ok
    real_t pde_mu_rs   = 1.0;                               // ok
    real_t pde_mu_cr   = pde_mu_r;                              // MISSING... Perhaps mu_r mispelled?
    real_t pde_k0      = 2*M_PI*pde_f*sqrt(pde_eps0*pde_mu0);  // ok

    std::complex<real_t> pde_eps_rs = pde_eps_rs_real + 1i*pde_eps_rs_imag;
    std::complex<real_t> pde_eps_crs = pde_eps_crs_real + 1i*pde_eps_crs_imag;

    std::complex<real_t> pde_bnd_const = 1.0/pde_eps_cr*(1i*pde_k0 + 1.0/(2*pde_r_t));

    std::complex<real_t> pde_eps_crs_inv = 1.0/pde_eps_crs;
    real_t pde_eps_cr_inv = 1.0/pde_eps_cr;

    gsMultiPatch<> pde_eps_cr_fun_real;
    gsMultiPatch<> pde_eps_cr_fun_imag;

    gsMultiPatch<> pde_mu_r_fun_real;
};




#endif //GSSTATEEQUATIONANTENNA_H
