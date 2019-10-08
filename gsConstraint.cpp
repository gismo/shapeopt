#include <gismo.h>
#include "gsConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

gsConstraint::gsConstraint(gsMultiPatch<>* mpin): m_mp(mpin){} // Only loads m_mp

void gsConstraint::evalCon_into(gsAsVector<real_t> & result)
{
    result = evalCon();
}

index_t gsConstraint::getSignOfPatch(index_t patch)
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsMultiPatch<> singlePatch(m_mp->patch(patch));
    gsMultiBasis<> dbasis(singlePatch);
    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(singlePatch);

    real_t integral = ev.integral(jac(G).det());

    return (integral > 0) - (integral < 0); // returning sign of integral
}
