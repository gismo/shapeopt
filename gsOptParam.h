#ifndef GSOPTPARAM_H
#define GSOPTPARAM_H
using namespace gismo;

#include "gsShapeOptProblem.h"
#include "gsSpringMethod.h"
#include "gsModLiao.h"
#include "gsWinslow.h"
#include "gsWinslowWithDeriv.h"
#include "gsLiao.h"
#include "gsHarmonic.h"
#include "gsOptParamMethod.h"
#include "gsAffineOptParamMethod.h"

class gsOptParam: public gsShapeOptProblem {
public:

    // Constructs from a multipatch, calling setupMappers, and setupOptParameters
    // Uses gsDetJacConstraint as default
    gsOptParam(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsMultiPatch<>> mp_goal, memory::shared_ptr<gsShapeOptLog> slog, index_t param, bool use_Lagrangian = false);

    gsOptParam(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsMultiPatch<>> mp_goal, std::vector< gsDofMapper > mappers, memory::shared_ptr<gsShapeOptLog> slog, index_t param, bool use_Lagrangian = false);

	void setup(index_t param, bool use_Lagrangian);

    // Evaluation of the objective, using the design contained in m_mp
    real_t evalObj() const ;

    // Evaluation of the gradient wrt design variables (tagged cps)
    gsVector<> gradObj() const ;

    // Evaluation of the gradient wrt all cps, needed for gsShapeOptWithReg
    gsVector<> gradAll() const ;

    void setupMappers();

    gsDofMapper mapper_grad() const;

    // Method to set the bounds on the design variables,
    void setupDesignBounds();

	void setQuad(real_t quA, index_t quB){ m_quA = quA; m_quB = quB; };

public:

	typedef memory::shared_ptr<gsOptParam> Ptr;

public:
    gsMultiPatch<>::Ptr m_mp_goal;
    gsSpringMethod m_pM_goal;
    gsVector<> m_tagged_goal;

	real_t m_quA = 4;
	real_t m_quB = 4;
};



# endif //GSOPTPARAM_H
