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
    gsOptParam(gsMultiPatch<>* mp, gsMultiPatch<>* mp_goal, gsShapeOptLog* slog, index_t param, bool use_Lagrangian = false);

    // Evaluation of the objective, using the design contained in m_mp
    real_t evalObj() const ;

    // Evaluation of the gradient wrt design variables (tagged cps)
    gsVector<> gradObj() const ;

    // Evaluation of the gradient wrt all cps, needed for gsShapeOptWithReg
    gsVector<> gradAll() const ;

    void setupMappers();

    gsDofMapper mapper_grad() const;

public:
    gsMultiPatch<> *m_mp_goal;
    gsSpringMethod m_pM_goal;
    gsVector<> m_tagged_goal;
};



# endif //GSOPTPARAM_H
