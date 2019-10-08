#include <gismo.h>
#include "gsConstraint.h"
#include "gsAggregatedConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

gsAggregatedConstraint::gsAggregatedConstraint(gsMultiPatch<>* mp, gsConstraint* constraint, gsFunction<>* map):
    gsConstraint(mp),
    m_mp(mp),
    m_constraint(constraint),
    m_map(map)
{
    setup();
}

gsVector<> gsAggregatedConstraint::evalCon()
{
    return m_map->eval(m_constraint->evalCon());
}


// FIXIT make this method faster!
gsIpOptSparseMatrix gsAggregatedConstraint::getJacobian()
{
    gsMatrix<> result;
    m_map->jacobian_into(m_constraint->evalCon(), result);
    gsIpOptSparseMatrix jac = m_constraint->getJacobian(); // Avoid conversion back and forth from gsIpOptSparseMatrix, preferably work directly with sparse matrices
    gsMatrix<> mat = result*jac.asDense(); // Avoid asDense() call
    gsIpOptSparseMatrix out(mat);
    return out;
}

gsVector<> gsAggregatedConstraint::getUpperBounds()
{
    gsVector<> out;
    out.setConstant(n_constraints, 1e19);
    return out;
}

gsVector<> gsAggregatedConstraint::getLowerBounds()
{
    gsVector<> out;
    out.setConstant(n_constraints, m_eps);
    return out;
}

// Implement
void gsAggregatedConstraint::setup()
{
    m_constraint->setup(); // Make sure the constraint is setup

    // Get m_space_mapper from m_constraint
    m_space_mapper = m_constraint->space_mapper();

    GISMO_ASSERT(m_constraint->n_constraints == m_map->domainDim(),"Constraint should have map.domainDim() constraints in gsAggregatedConstraint..");

    // Number of constraints is equal the target dimension of the map
    n_constraints = m_map->targetDim();
}
