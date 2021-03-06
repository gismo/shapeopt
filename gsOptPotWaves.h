#ifndef GSOPTPOTWAVES_H
#define GSOPTPOTWAVES_H
using namespace gismo;

#include "gsShapeOptProblem.h"
#include "gsStateEquationPotWaves.h"
#include "gsSpringMethod.h"
#include "gsModLiao.h"
#include "gsWinslow.h"
#include "gsWinslowWithDeriv.h"
#include "gsLiao.h"
#include "gsHarmonic.h"
#include "gsOptParamMethod.h"
#include "gsAffineOptParamMethod.h"

class gsOptPotWaves: public gsShapeOptProblem {
public:

    // Constructs from a multipatch, calling setupMappers, and setupOptParameters
    // Uses gsDetJacConstraint as default
    gsOptPotWaves(memory::shared_ptr<gsMultiPatch<>> mp, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, index_t param, real_t quA, index_t quB, bool useDetJCons = false, bool use_Lagrangian = false);

    // FIXIT implemtn a way to have mp_ptr defined by this class
    //gsOptPotWaves(index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, index_t param, real_t quA, index_t quB, bool useDetJCons = false, bool use_Lagrangian = false);

    void constructor(index_t param, real_t quA, index_t quB, bool use_Lagrangian);

    // Evaluation of the objective, using the design contained in m_mp
    real_t evalObj() const ;

    // Evaluation of the gradient wrt design variables (tagged cps)
    // Main implementations are in gradAll()
    gsVector<> gradObj() const ;

    // Evaluation of the gradient wrt all control points
    gsVector<> gradAll() const ;

    // Solve if u != m_state
    void solveStateEquation( ) const ;

    void getObjVec() const;
    void getGradObjMat() const;

    // Mappers for the problem,
    // Interfaces are glued, the pml patches is eliminated, Gamma_s is tagged and eliminated.
    void setupMappers();

    void setupMapperGrad();

    // Method to set the bounds on the design variables,
    void setupDesignBounds();

    void setupGoalFunctions();
    void plotGoalFunctions(std::string name);

    gsVector<> getObjDerivativeDc(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;
    gsVector<> getObjDerivativeDu(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;

    gsDofMapper mapper_grad() const;

    bool isFlatInPML(index_t i);
    bool isPatchInDomain(index_t p);
    bool isCpsInDomain(index_t i, index_t p, index_t dim);

    void eliminatePatch(index_t p);

    void getTestMatrix();

    gsMultiPatch<> getUR() { return m_uR; };
    gsMultiPatch<> getUI() { return m_uI; };


public:
    // Define smart pointers
	typedef memory::shared_ptr<gsOptPotWaves> Ptr;
	typedef memory::unique_ptr<gsOptPotWaves> uPtr;

public:
    // State equation class to hold state equation
    mutable gsStateEquationPotWaves m_stateEq;

    real_t desLowerBoundx = -0.75;
    real_t desUpperBoundx = 0.75;

    real_t desLowerBoundy = -0.5;
    real_t desUpperBoundy = 0.75;

    real_t desLowerBoundz = -0.5;

    // Mapper used for assembling gradients
    mutable gsDofMapper m_mapper_grad;
    mutable bool m_mapper_grad_exist = false;

    index_t m_dim = 3;
    index_t n_controlpoints;

    mutable gsVector<> m_objVec_re;
    mutable gsVector<> m_objVec_im;

    mutable gsMatrix<> m_gradObjMat_re;
    mutable gsMatrix<> m_gradObjMat_im;

    gsMatrix<> m_testMatrix;

    // Hold latest solve of state equation
    mutable gsMultiPatch<> m_uR, m_uI;
    mutable gsVector<> m_uRVec, m_uIVec;
    mutable gsVector<> m_state; // Save the design used to calculate u_real and u_imag

    gsFunctionExpr<>::uPtr goal_fun_re;
    gsFunctionExpr<>::uPtr goal_fun_im;
    gsFunctionExpr<>::uPtr weight_fun;
                  
    gsFunctionExpr<>::uPtr grad_goal_fun_re;
    gsFunctionExpr<>::uPtr grad_goal_fun_im;
    gsFunctionExpr<>::uPtr grad_weight_fun;

};



# endif //GSOPTPOTWAVES_H
