#ifndef GSOPTANTENNA_H
#define GSOPTANTENNA_H
using namespace gismo;

#include "gsShapeOptProblem.h"
#include "gsStateEquationAntenna.h"
#include "gsSpringMethod.h"
#include "gsModLiao.h"
#include "gsWinslow.h"
#include "gsWinslowWithDeriv.h"
#include "gsLiao.h"
#include "gsHarmonic.h"
#include "gsOptParamMethod.h"
#include "gsAffineOptParamMethod.h"

class gsOptAntenna: public gsShapeOptProblem {
public:

    // Constructs from a multipatch, calling setupMappers, and setupOptParameters
    // Uses gsDetJacConstraint as default
    gsOptAntenna(memory::shared_ptr<gsMultiPatch<>> mp, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, index_t param, real_t quA, index_t quB, bool useDetJCons = false, bool use_Lagrangian = false);

    // Constructs from a multipatch, calling setupMappers, and setupOptParameters
    // Uses specified constraint from \a constraint
    gsOptAntenna(memory::shared_ptr<gsMultiPatch<>> mp, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, memory::shared_ptr<gsConstraint> constraint, index_t param, real_t quA, index_t quB, bool use_Lagrangian = false);

    // Evaluation of the objective, using the design contained in m_mp
    real_t evalObj() const ;

    // Evaluation of the gradient wrt design variables (tagged cps)
    // Main implementations are in gradAll()
    gsVector<> gradObj() const ;

    // Evaluation of the gradient wrt all control points
    gsVector<> gradAll() const ;

    // Mappers for antenna problem,
    // Interfaces are glued, the bnd of patch 3 is tagged and eliminated.
    void setupMappers();

    // Method to set the bounds on the design variables,
    void setupDesignBounds();

    // Helper methods to calculate the gradient
    gsVector<> evaluateDerivativeTerm1(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;
    gsVector<> evaluateDerivativeTerm2(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;
    gsVector<> gradientObjWithoutAdjoint() const;
    gsVector<> getObjDerivativeDu(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;

    gsDofMapper mapper_grad() const;

    real_t getQuA(){ return m_stateEq.m_quA; };
    real_t getQuB(){ return m_stateEq.m_quB; };
public:

	typedef memory::shared_ptr<gsOptAntenna> Ptr;

public:
    // function for objective function and its derivatives
    gsFunctionExpr<>::Ptr delta;
    gsFunctionExpr<>::Ptr ddeltadx;
    gsFunctionExpr<>::Ptr ddeltady;

    // State equation class to hold state equation
    mutable gsStateEquationAntenna m_stateEq;

    // Design bounds
    real_t d_g = 57.5*1e-9;       // Minimum allowed gap between the two antennas 57.5 [nm]
    real_t d_bw = 776.25*1e-9;    // Width of bounding box (Around x=0)
    real_t d_bh = 546.25*1e-9;    // Height of bounding box

    // Mapper used for assembling gradients
    mutable gsDofMapper m_mapper_grad;
    mutable bool m_mapper_grad_exist = false;

    // The patch of the antenna
    index_t fixedPatch = 3;


};



# endif //GSOPTANTENNA_H
