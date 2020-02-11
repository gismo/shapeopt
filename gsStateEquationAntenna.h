#ifndef GSSTATEEQUATIONANTENNA_H
#define GSSTATEEQUATIONANTENNA_H

#include <math.h>
#include <complex>
#include <sstream>
#include <string>
using namespace gismo;

// FIXIT: Clean up. Make more flexible wrt to topology, e.g. use Piecewise function
class gsStateEquationAntenna{
public:
    gsStateEquationAntenna(gsMultiPatch<>* mpin, index_t numRefine): mp(mpin), dbasis(*mp), zero("0.0",2){
        dbasis.setDegree(degree);
        // gsInfo << "OBS CHECK DISCRETIZATION OF PDE!!\n";
        for (index_t i = 0; i < numRefine; i++){
            dbasis.uniformRefine();
        }

	std::complex<real_t> _i(0,1);

    	pde_eps_rs = pde_eps_rs_real + _i*pde_eps_rs_imag;
    	pde_eps_crs = pde_eps_crs_real + _i*pde_eps_crs_imag;

    	pde_bnd_const = 1.0/pde_eps_cr*(_i*pde_k0 + 1.0/(2*pde_r_t));

        pde_f       = pde_f/pde_L_f;
        pde_omega   = 2*M_PI*pde_f;
        pde_eps_rs_real  = -20.198689873114;
        pde_eps_rs_imag = 1.38170639237;
        pde_eps_crs_real = pde_eps_rs_real;
        pde_eps_crs_imag = pde_eps_rs_imag - pde_sigma/(pde_omega*pde_L_f*pde_eps0);
        pde_k0      = 2*M_PI*pde_f*sqrt(pde_eps0*pde_mu0);

        pde_eps_rs = pde_eps_rs_real + _i*pde_eps_rs_imag;
        pde_eps_crs = pde_eps_crs_real + _i*pde_eps_crs_imag;

        pde_bnd_const = 1.0/pde_eps_cr*(_i*pde_k0 + 1.0/(2*pde_r_t));

        pde_eps_crs_inv = 1.0/pde_eps_crs;
        pde_eps_cr_inv = 1.0/pde_eps_cr;

        //FIXIT Antenna patch is hardcoded to 3.!!!
        pde_eps_cr_fun_real = getPieceWiseFunctionOnPatch(mp->nBoxes(),3,pde_eps_crs_inv.real(), pde_eps_cr_inv);
        pde_eps_cr_fun_imag = getPieceWiseFunctionOnPatch(mp->nBoxes(),3,pde_eps_crs_inv.imag(), 0.0);
        // pde_eps_cr_fun_real = getPieceWiseFunctionOnPatch(mp->nBoxes(),3,pde_eps_cr_inv, pde_eps_cr_inv);
        // pde_eps_cr_fun_imag = getPieceWiseFunctionOnPatch(mp->nBoxes(),3,0.0, 0.0);
        pde_mu_r_fun_real = getPieceWiseFunctionOnPatch(mp->nBoxes(),3,pde_mu_rs, pde_mu_r);
        // plotSolution(pde_eps_cr_fun_real,"pde_eps_cr_fun_real");
        // plotSolution(pde_eps_cr_fun_imag,"pde_eps_cr_fun_imag");

        // printConstants();

        // Change domain size
        // for (index_t i = 0; i < mp->nPatches(); i++){
        //    gsMatrix<> cc = mp->patch(i).coefs();
        //    mp->patch(i).setCoefs(pde_r_t/8.0*cc);
        // }

        char tmp[200];
        snprintf(tmp, 200, "(-%f*sqrt(%f*%f))", pde_k0, pde_eps_cr, pde_mu_r);
        std::string k_str = tmp;

        snprintf(tmp, 200, "cos(%s*x)", k_str.c_str());
        std::string Hiz_real_str = tmp;
        // gsInfo << "\n Hiz real: " << Hiz_real_str << "\n";
        Hiz_real = *(new gsFunctionExpr<>(Hiz_real_str,2));

        snprintf(tmp, 200, "sin(%s*x)", k_str.c_str());
        std::string Hiz_imag_str = tmp;
        // gsInfo << "\n Hiz imag: " << Hiz_imag_str << "\n";
        Hiz_imag = *(new gsFunctionExpr<>(Hiz_imag_str,2));

        snprintf(tmp, 200, "(-%s*sin(%s*x))", k_str.c_str(), k_str.c_str());
        std::string dHizdx_real_str = tmp;
        // gsInfo << "\n dHizdx real: " << dHizdx_real_str << "\n";
        dHizdx_real = *(new gsFunctionExpr<>(dHizdx_real_str,2));

        snprintf(tmp, 200, "(%s*cos(%s*x))", k_str.c_str(), k_str.c_str());
        std::string dHizdx_imag_str = tmp;
        // gsInfo << "\n dHizdx imag: " << dHizdx_imag_str << "\n";
        dHizdx_imag = *(new gsFunctionExpr<>(dHizdx_imag_str,2));

        snprintf(tmp, 200, "x/%f*%s", pde_r_t, dHizdx_real_str.c_str());
        std::string dHizdn_real_str = tmp;
        // gsInfo << "\n dHizdn real: " << dHizdn_real_str << "\n";
        dHizdn_real = *(new gsFunctionExpr<>(dHizdn_real_str,2));

        snprintf(tmp, 200, "x/%f*%s", pde_r_t, dHizdx_imag_str.c_str());
        std::string dHizdn_imag_str = tmp;
        // gsInfo << "\n dHizdn imag: " << dHizdn_imag_str << "\n";
        dHizdn_imag = *(new gsFunctionExpr<>(dHizdn_imag_str,2));

        snprintf(tmp, 200, "1/%f*%s - x/%f*%s^2*%s", pde_r_t, dHizdx_real_str.c_str(), pde_r_t, k_str.c_str(), Hiz_real_str.c_str());
        std::string d2Hizdndx_real_str = tmp;
        // gsInfo << "\n d2Hizdndx real: " << d2Hizdndx_real_str << "\n";
        d2Hizdndx_real = *(new gsFunctionExpr<>(d2Hizdndx_real_str,2));

        snprintf(tmp, 200, "1/%f*%s - x/%f*%s^2*%s", pde_r_t, dHizdx_imag_str.c_str(), pde_r_t, k_str.c_str(), Hiz_imag_str.c_str());
        std::string d2Hizdndx_imag_str = tmp;
        // gsInfo << "\n dHizdn imag: " << d2Hizdndx_imag_str << "\n";
        d2Hizdndx_imag = *(new gsFunctionExpr<>(d2Hizdndx_imag_str,2));

        bcInfo.addCondition(4, boundary::south, condition_type::neumann, &zero);
        // bcInfoZero.addCondition(4, boundary::east, condition_type::neumann, &zero); //zero bcs
        // bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &Hiz_real);
        // bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, &Hiz_real);
        // bcInfo.addCondition(2, boundary::west,  condition_type::dirichlet, &Hiz_real);
        // bcInfo.addCondition(4, boundary::east, condition_type::dirichlet, &Hiz_real);
    };

    gsMatrix<> getU(index_t realOrImag);
    gsMatrix<> getAu(index_t realOrImag, gsMatrix<> U);

    gsMatrix<> getDerivativeOfAu(index_t realOrImag, gsMultiPatch<> sol){
        return getDerivativeOfAuPart1(realOrImag,sol) + getDerivativeOfAuPart2(realOrImag,sol);
    };
    gsMatrix<> getDerivativeOfAuPart1(index_t realOrImag, gsMultiPatch<> sol); // For implenting new gradients
    gsMatrix<> getDerivativeOfAuPart2(index_t realOrImag, gsMultiPatch<> sol);


    gsMatrix<> getDerivativeWithoutSolving(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag);
    gsVector<> solveAdjoint(gsVector<> &rhs);

    gsMatrix<> solve(real_t &err, index_t &numDofs, index_t &numElems, index_t numRefine);
    gsMatrix<> getDerivativeOfRhsZeroBC(index_t realOrImag);
    gsMatrix<> getRhsZeroBC(index_t realOrImag);
    gsMatrix<> getDerivativeOfKu(gsMultiPatch<> sol);
    gsMatrix<> getKu(gsMultiPatch<> sol);
    gsMatrix<> getDerivativeOfU();
    gsMultiPatch<> solve();
    gsMatrix<> getU();
    real_t getResidual();
    void plotMesh(gsMatrix<> solVector);
    void printMatSize(gsMatrix<> mat, std::string name);
    void getFandMS(gsFunctionExpr<> *&f, gsFunctionExpr<> *&ms, gsFunctionExpr<> *&df_dx, gsFunctionExpr<> *&df_dy);
    void getMSDerivatives(gsFunctionExpr<> *&dms_dx, gsFunctionExpr<> *&dms_dy);
    gsMultiPatch<> getPieceWiseFunctionOnPatch(index_t nBoxes, index_t patch, real_t val_on_patch, real_t val_elsewhere);

    void getTerm(index_t realOrImag, gsSparseMatrix<> &mat, gsVector<> &rhs);
    void getSystem(gsSparseMatrix<> &A, gsVector<> &rhs);
    void assembleAndSolve();
    void solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag);
    void plotSolution(gsMultiPatch<> &sol, std::string name);
    void plotSolution(std::string name);
    void plotMagnitude(std::string name);

    void printConstants();

    void writeToFile(gsMatrix<> mat, std::string name) const{
        gsInfo << "WRITING to " << name << "\n";
        std::ofstream f(name);
        for(index_t i = 0; i < mat.rows(); i++){
            for(index_t j = 0; j < mat.cols(); j++){
                f << std::setprecision(20) << mat(i,j) << " ";
            }
            f << "\n";
        }
    }

    real_t quA() { return m_quA; };
    index_t quB() { return m_quB; };


public:
    gsMultiPatch<>* mp;
    index_t degree = 2;
    bool isRefined = false;
    index_t numRef = 1;

    // Boundary conditions
    gsFunctionExpr<> zero;

    gsFunctionExpr<> Hiz_real;
    gsFunctionExpr<> Hiz_imag;

    gsFunctionExpr<> dHizdx_real;
    gsFunctionExpr<> dHizdx_imag;

    gsFunctionExpr<> dHizdn_real;
    gsFunctionExpr<> dHizdn_imag;
    gsFunctionExpr<> d2Hizdndx_real;
    gsFunctionExpr<> d2Hizdndx_imag;

    gsBoundaryConditions<> bcInfo;
    gsBoundaryConditions<> bcInfoZero;

    // Assembler variables
    gsMultiBasis<> dbasis;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsSparseSolver<>::LU solver;
    gsMatrix<> solVector;

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

    std::complex<real_t> pde_eps_rs;
    std::complex<real_t> pde_eps_crs;

    std::complex<real_t> pde_bnd_const;

    std::complex<real_t> pde_eps_crs_inv = 1.0/pde_eps_crs;
    real_t pde_eps_cr_inv = 1.0/pde_eps_cr;

    gsMultiPatch<> pde_eps_cr_fun_real;
    gsMultiPatch<> pde_eps_cr_fun_imag;

    gsMultiPatch<> pde_mu_r_fun_real;

    real_t m_quA = 2;   // Quadrature option
    index_t m_quB = 2;   // Quadrature option
};




#endif //GSSTATEEQUATIONANTENNA_H
