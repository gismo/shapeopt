// FIXIT: Write intro..
#ifndef GSAFFINEOPTPARAMMETHOD_H
#define GSAFFINEOPTPARAMMETHOD_H
using namespace gismo;

#include "gsAffineParamMethod.h"

// FIXIT: Think of a better name for this class? QuadProg?
class gsAffineOptParamMethod : public gsAffineParamMethod{
public:
    // Construct from gsOptParamMethod
    gsAffineOptParamMethod(gsOptParamMethod* optParamMethod);

    // Method to get update as a vector, to be used to calculate A and b
    //      from gsAffineParamMethod
    // Returns: vector of free control points
    gsVector<> getUpdate(gsVector<> x);

    // Evaluation of affine objective function
    real_t evalObj(real_t obj, gsVector<> grad, gsMatrix<> hess, gsVector<> x, gsVector<> x0);

    // Updates the reference parametrization to the one currently hold in m_problem;
    void reset();

    // Accesors
    index_t numDesignVars(){ return n_free; }
    gsVector<> refFree(){ return m_refFree; }
    gsVector<> refTagged(){ return m_refTagged; }

protected:
    // Pointer to gsOptParamMethod
    gsOptParamMethod* m_optParamMethod;

    // Control points from reference parametrization
    gsVector<> m_refFree;
    gsVector<> m_refTagged;

    // Entities that define the map
    // grad and hess with respect to free variables
    gsVector<> m_grad;
    gsMatrix<> m_hess;
    real_t obj;

    // FIXIT: Is there a sparse structure to exploit here?
    // KKT system
    gsMatrix<> m_KKTsystem;
    gsVector<> m_rhs;

    // Solver to solve KKT system
    Eigen::FullPivLU<Eigen::Matrix<real_t,Dynamic,Dynamic>> m_solver;
};


#endif //GSAFFINEOPTPARAMMETHOD_H
