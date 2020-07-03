#ifndef GSSTATEEQUATION_H
#define GSSTATEEQUATION_H

#include <math.h>
#include <complex>
#include <sstream>
#include <string>
using namespace gismo;

class gsStateEquation{
public:
    gsStateEquation(gsMultiPatch<>::Ptr mpin, index_t numRefine): m_mp(mpin), dbasis(*m_mp){
        dbasis.setDegree(degree);
        // gsInfo << "OBS CHECK DISCRETIZATION OF PDE!!\n";
        for (index_t i = 0; i < numRefine; i++){
            dbasis.uniformRefine();
        }
    };

    gsStateEquation() 
    { 
        // Empty constructor
    };

    // Solves Au = F, and returns either real or imaginary part
    gsMatrix<> getU(index_t realOrImag);

    // Method to multiply A with arbitrary matrix u.
    gsMatrix<> getAu(index_t realOrImag, gsMatrix<> U);

    // 
    // Methods for calculation of the derivative, generic methods (not PDE specific)
    //
    
    // Method to get du/dc 
    //
    // u is given by the linear system of equation
    // Au = F
    //
    // Its derivative is calculated as
    // du/dc = A^{-1} ( dF/dc - dA/dc*u )
    gsMatrix<> getDerivativeOfU();

    // Method to get A*du/dc 
    // A*du/dc =  dF/dc - dA/dc*u 
    // so the methods simply returns dF/dc - dA/dc*u
    // can be used for the adjoint method
    gsMatrix<> getDerivativeWithoutSolving(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag);

    // Solves the linear system
    // Au = rhs
    // with arbitrary rhs.
    // Can be used in an adjoint method
    gsVector<> solveAdjoint(gsVector<> &rhs);

    // Method to solve PDE and return real and imaginary part in u_real and u_imag
    virtual void solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag);
    virtual void solve(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag, gsVector<> &solVec_real, gsVector<> &solVec_imag);

    // Method to get real and imaginary part of F
    gsMatrix<> getRhsZeroBC(index_t realOrImag);

    // Method to get solution as vector (both real and imaginary part)
    gsMatrix<> getU();

    // Get the full real system build from the imaginary and real parts.
    // Note that the system is symmetric
    // [A_r -A_i
    //  -A_i -A_r]
    void getSystem(gsSparseMatrix<> &A, gsVector<> &rhs);

    //
    // To be overloaded
    // 

    // Method to get real and imaginary part of dA/dc*u
    // Takes 'u' as a gsMultiPatch as an input
    virtual gsMatrix<> getDerivativeOfAu(index_t realOrImag, gsMultiPatch<> sol)
    { GISMO_NO_IMPLEMENTATION; };

    // Method to get real and imaginary part of  dF/dc
    virtual gsMatrix<> getDerivativeOfRhsZeroBC(index_t realOrImag) { GISMO_NO_IMPLEMENTATION; };

    // Method to get real and imaginary part of the system. The result is stored in mat, and rhs
    // Note that this method is pure virtual so it has to be overloaded
    virtual void getTerm(index_t realOrImag, gsSparseMatrix<> &mat, gsVector<> &rhs) = 0;

    virtual void printConstants() { GISMO_NO_IMPLEMENTATION; };

    // 
    // Methods to visualize content and solutions
    //

    void printMatSize(gsMatrix<> mat, std::string name);
    void plotSolution(gsMultiPatch<> &sol, std::string name);
    void plotSolution(std::string name);
    void plotMagnitude(std::string name);

    void writeToFile(gsMatrix<> mat, std::string name) const;

    // 
    // Helper methods
    //

    gsMultiPatch<> getPieceWiseFunctionOnPatch(index_t nBoxes, index_t patch, real_t val_on_patch, real_t val_elsewhere);



    // 
    // Accessors
    //
    real_t quA() { return m_quA; };
    index_t quB() { return m_quB; };

    void setQuad(real_t A, index_t B){ m_quA = A; m_quB = B; };


public:
    // Holder for last solution
    gsMultiPatch<> m_ur; 
    gsMultiPatch<> m_ui; 

    // Physical domain
    gsMultiPatch<>::Ptr m_mp;

    // Degree of splines used for analysis
    index_t degree = 2;

    // How many times should the 'geometry' splines be refined for analysis. 
    index_t m_numRefine = 1;
    bool isRefined = false;

    // Assembler variables
    gsMultiBasis<> dbasis;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsSparseSolver<>::LU solver;
    gsMatrix<> solVector;

    real_t m_quA = 2;   // Quadrature option
    index_t m_quB = 2;   // Quadrature option
};




#endif //GSSTATEEQUATION_H
