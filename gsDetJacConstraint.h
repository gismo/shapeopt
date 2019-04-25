#ifndef DETJACCONSTRAINT_H
#define DETJACCONSTRAINT_H

#include "IpOptSparseMatrix.h"
using namespace gismo;

class liaoOptProblem;

class detJacConstraint{
public:
  detJacConstraint(gsMultiPatch<>* mpin);

	gsVector<> generateDResultVector();

	void getDvectors(gsVector<> &result);

  gsVector<> getDvectors();

  void getDerivRhsFromPatch(index_t patch, gsSparseMatrix<> &xJac, gsSparseMatrix<> &yJac);

  void getJacobianFromPatch(index_t patch, gsMatrix<> &xJac, gsMatrix<> &yJac);

  IpOptSparseMatrix getJacobian();

  gsVector<> getDesignVariables();

  void updateDesignVariables(gsVector<> des);

  const gsMultiBasis<> & detJacBasis() const { return m_detJacBasis; }
  const bool & isSolverSetup(index_t i) const { return m_areSolversSetup[i]; }
  friend class liaoOptProblem;

  gsVector<> getUpperBounds(real_t eps);

  void plotDetJ(std::string name);
public:
    gsMultiPatch<>* mp;
	gsMultiBasis<> m_detJacBasis;

	gsVector<gsSparseSolver<>::LU> solversMassMatrix;
	gsVector<bool> m_areSolversSetup;

    index_t n_controlpoints;
    index_t n_constraints;
};




#endif //DETJACCONSTRAINT_H
