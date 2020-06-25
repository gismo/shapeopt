#ifndef GSSTATEEQUATIONPOTWAVES_H
#define GSSTATEEQUATIONPOTWAVES_H

#include <math.h>
#include <complex> 
#include <sstream> 
#include <string>
#include "gsStateEquation.h"

using namespace gismo;

class gsStateEquationPotWaves: public gsStateEquation {
public:
    gsStateEquationPotWaves(index_t numRefine);

    gsStateEquationPotWaves(gsMultiPatch<>::Ptr mp_ptr, index_t numRefine);

    gsStateEquationPotWaves(index_t numRefine, real_t Lx, real_t Ly, real_t Lz, real_t lx, real_t ly, real_t lz); 

    void constructor();

    // Method to get the real and imaginary part of the linear system and rhs
    // Real part : realOrImag == 0
    // Imag part : realOrImag == 1
    // 
    // Result is saved in arguments mat and rhs.
    void getTerm(index_t realOrImag, gsSparseMatrix<> &mat, gsVector<> &rhs);

    // Overloaded to allow the use of Dirichlet bnd conditions
    void solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag);

    // Method to mark boundaries Gamma_f and Gamma_s
    // Result go in members bcInfo_Gamma_f and bcInfo_Gamma_s
    void markBoundaries();
    void markBoundariesDirichlet();
    void markBoundariesDirichletNoPML();

    void markBoundariesDirichletNoPML_NoCenter();

    // Method to generate initial domain.
    gsMultiPatch<>::Ptr getInitialDomain(bool includeCenter = false);
    gsMultiPatch<>::Ptr getInitialDomainTmp(bool includeCenter = false);
    gsMultiPatch<>::Ptr getInitialDomainNoPML(bool includeCenter = false);
    gsMultiPatch<>::Ptr getInitialDomainNoPMLInZDir(bool includeCenter = false);

    void getUI(gsMultiPatch<> &uI_re, gsMultiPatch<> &uI_im) ;

    void plotVelocityField(gsMultiPatch<> &ur, gsMultiPatch<> &ui, real_t timestep, std::string outfile, bool includeIncident = true);

    // Convergence test against manufactures solution
    void convergenceTest(index_t maxRefine, std::string outfolder); 

    // Convergence test against manufactures solution
    void convergenceTestNoPML_NoCenter(index_t maxRefine, std::string outfolder); 

    // Convergence test against manufactures solution, center filled with water (to only test PML)
    void convergenceTestOnlyPML(index_t maxRefine, std::string outfolder); 

    // Convergence test against manufactures solution, center filled with water (to only test PML)
    void convergenceTestNoPML(index_t maxRefine, std::string outfolder); 

    // Convergence test against manufactures solution, center filled with water (to only test PML)
    void convergenceTestOnlyPMLAllDir(index_t maxRefine, std::string outfolder, real_t angle = 0); 

    void convTestHelper(bool useDirichlet, bool useNeumann, index_t max_refine, std::string outfolder);

    void pointSourceTest(std::string outfolder);

    void pointSourceTestForce(std::string outfolder);


    void setup();

    gsMultiPatch<> getPieceWiseFunctionOnPatch(gsVector<> vals);

    gsMultiPatch<> markPatches(gsVector< index_t > patches);

    gsMultiPatch<> markInnerDomain();

    void testSplineSpace(real_t k);

    bool isBndGamma_s( index_t b); // Input is the boundary index!

    // Methods for derivatives
    //
    gsMatrix<> getDerivativeOfRhsZeroBC(index_t realOrImag);
    gsMatrix<> getDerivativeOfAu(index_t realOrImag, gsMultiPatch<> sol);

    bool isPatchInDomain(index_t p);

    gsVector< index_t > getVectorWithDomainPatches();
    gsVector< index_t > getVectorWithPMLPatches();


private:
    // Helper method
    gsMatrix<> getCoefs(index_t p);

    index_t m_numRefine;


public:

    //gsSparseSolver<>::CGDiagonal solver;
    gsMultiBasis<> dbasis_domain;
    gsMultiPatch<> m_mp_domain;

    // Wave parameters
    real_t wave_omega;//= sqrt(K * g) set in constructor
    real_t wave_g       = 9.82;                         // Gravitation
    real_t wave_K       = 4.0;   // Wave number
    real_t wave_A       = 1.0;                            // Amplitude

    // PML parameters
    real_t pml_n        = 3;    // Exponent
    real_t pml_C        = 5;    // constant

    real_t pml_Lx       = 4;
    real_t pml_lx       = 3;
    real_t pml_Ly       = 4;
    real_t pml_ly       = 3;
    real_t pml_Lz       = 3;
    real_t pml_lz       = 2.0;

    real_t init_lbx     = 0.5;
    real_t init_lby     = 0.5;
    real_t init_lbz     = 0.5;
    //real_t pml_Lx       = 8;
    //real_t pml_lx       = 1;
    //real_t pml_Ly       = 8;
    //real_t pml_ly       = 2*M_PI;
    //real_t pml_Lz       = 1.2;
    //real_t pml_lz       = 1;

    //real_t init_lbx     = 0.5;
    //real_t init_lby     = 2.0;
    //real_t init_lbz     = 0.5;

    std::map< std::string, real_t > const_map;

    gsFunctionExpr<>::uPtr wave_u_I_re;
    gsFunctionExpr<>::uPtr wave_u_I_im;

    gsFunctionExpr<>::uPtr strech_Jm1_detJinv_asVec_re;
    gsFunctionExpr<>::uPtr strech_Jm1_detJinv_asVec_im;

    gsFunctionExpr<>::uPtr strech_detJ_Jm1_Jm1_asVec_re;
    gsFunctionExpr<>::uPtr strech_detJ_Jm1_Jm1_asVec_im;

    gsFunctionExpr<>::uPtr strech_detJ_re;
    gsFunctionExpr<>::uPtr strech_detJ_im;

    gsFunctionExpr<>::uPtr strech_detJ_len_re;
    gsFunctionExpr<>::uPtr strech_detJ_len_im;

    gsFunctionExpr<> zero;

    // <F,n> is the rhs of Neumann condition at Famma_s
    // in the state equation it is -grad(u_I)
    // For convergence test with manufactured solutions
    // F = grad(u_sol), with u_sol = cos(Ky)exp(Kz)
    gsFunctionExpr<>::uPtr pde_F_re;
    gsFunctionExpr<>::uPtr pde_F_im;

    gsFunctionExpr<>::uPtr pde_dF_re;
    gsFunctionExpr<>::uPtr pde_dF_im;

    // Gaussian bell curve for testing with point source
    bool                   useNeuBell = false;
    real_t                 bell_sigma       = 0.1;
    real_t                 bell_amplitude   = 1;
    gsFunctionExpr<>::uPtr bell_term_re;
    gsFunctionExpr<>::uPtr bell_term_im;

    // Force term for testing with point source
    bool                   useForce = false;
    gsFunctionExpr<>::uPtr lapl_uEx;
                                      
    // Boundary conditions
    gsBoundaryConditions<> bcInfo_Gamma_f; // Marks free surface, Gamma_f, as neumann
    gsBoundaryConditions<> bcInfo_Gamma_b; // Marks bottom surface, Gamma_b, as neumann
    gsBoundaryConditions<> bcInfo_Gamma_s; // Marks boundary of cloacking device, Gamma_s, as neumann

    gsBoundaryConditions<> bcInfo_Dirichlet_re;   // Marks Dirichlet (y=0) bnd for convergence testing
    gsBoundaryConditions<> bcInfo_Dirichlet_im;   // Marks Dirichlet (y=0) bnd for convergence testing
    gsFunctionExpr<>::uPtr pde_u_dir_re;          // Solution at Dirichlet bnd
    gsFunctionExpr<>::uPtr pde_u_dir_im;          // Solution at Dirichlet bnd HAS TO BE ZERO

    bool useDir = false;
    bool useNeu = true;

    gsVector< index_t > m_isBndGamma_s_vec;
    
    index_t m_dim = 3;

};

#endif //GSSTATEEQUATIONPOTWAVES_H
