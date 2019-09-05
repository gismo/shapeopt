/** @file gsAffineOptParamMethod.h

@brief  Implements an quadratic program, that approximates an optimization problems
        from gsOptParamMethod. Basicly it can be seen as the linearization of a
        nonlinear optimization based parametrization strategy.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/

#ifndef GSAFFINEOPTPARAMMETHOD_H
#define GSAFFINEOPTPARAMMETHOD_H
using namespace gismo;

#include "gsOptParamMethod.h"
#include "gsAffineParamMethod.h"

// FIXIT: Think of a better name for this class? QuadProg?
class gsAffineOptParamMethod : public gsAffineParamMethod{
public:
    // Construct from gsOptParamMethod
    gsAffineOptParamMethod(gsOptParamMethod* optParamMethod, bool use_Lagrangian = false);

    // Method to get update as a vector, to be used to calculate A and b
    //      from gsAffineParamMethod
    // Returns: vector of free control points
    gsVector<> getUpdate(gsVector<> x);

    // Evaluation of affine objective function
    // c should hold the free cps and x the tagged
    real_t evalObj(gsVector<> c, gsVector<> x);

    // Update Reference parametrization by first solving m_optParamMethod and then reset
    void updateAndReset() { m_optParamMethod->update(); reset(); };

    // Updates the reference parametrization to the one currently hold in m_problem;
    // FIXIT:   This should perhaps compute also a new parametrization by calling
    //          gsOptParamMethod.update();
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
    real_t m_obj;
    gsVector<> m_grad;
    gsMatrix<> m_hess;
    gsMatrix<> m_hessTagged; // Holds the second order derivatives wrt free and tagged

    // FIXIT: Is there a sparse structure to exploit here?
    // KKT system
    gsMatrix<> m_KKTsystem;
    gsVector<> m_rhs;

    // Solver to solve KKT system
    Eigen::FullPivLU<Eigen::Matrix<real_t,Dynamic,Dynamic>> m_solver;

    // Boolean that is true if the linearization should be performed
    // on the Lagrangian instead of the objective of the original
    // optimization
    bool m_use_Lagrangian;

};


#endif //GSAFFINEOPTPARAMMETHOD_H
