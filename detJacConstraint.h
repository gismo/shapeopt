#ifndef DETJACCONSTRAINT_H
#define DETJACCONSTRAINT_H

#include "IpOptSparseMatrix.h"
using namespace gismo;

class liaoOptProblem;

class detJacConstraint{
public:
  detJacConstraint(gsMultiPatch<>* mpin);

	void getDvectors(gsVector<> &result);

	gsVector<> getDvectors();

  IpOptSparseMatrix getJacobian();

  gsVector<> getDesignVariables();

  void updateDesignVariables(gsVector<> des);

  const gsMultiBasis<> & detJacBasis() const { return m_detJacBasis; }
  const bool & isSolverSetup() const { return m_isSolverSetup; }
  friend class liaoOptProblem;

  gsVector<> getUpperBounds(real_t eps);

  void plotDetJ(std::string name);
public:
  gsMultiPatch<>* mp;
	gsMultiBasis<> m_detJacBasis;

	gsSparseSolver<>::LU solverMassMatrix;
	bool m_isSolverSetup;

  gsVector<> M,Mm1;

  index_t n_controlpoints;
  index_t n_constraints;
};




#endif //DETJACCONSTRAINT_H
