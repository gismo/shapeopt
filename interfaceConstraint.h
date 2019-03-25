#ifndef INTERFACECONSTRAINT_H
#define INTERFACECONSTRAINT_H
using namespace gismo;

class liaoOptProblem;

#include "IpOptSparseMatrix.h"

class interfaceConstraint{
public:

  interfaceConstraint(gsMultiPatch<>* mpin);

  gsMatrix<> generateConstraintMatrix();

  void addInterfaceConstraint(index_t i, index_t j, index_t &count, gsMatrix<> &interfaceConstraintMatrix);

  gsVector<index_t> getPatchSideIndicies(patchSide ps);

  void getLocalPatchSideIndicies(patchSide ps, gsVector<index_t> &patchIndicies);

  gsVector<> getUpperBounds();

  gsVector<> getUpperBounds(index_t fixedPatch);

  gsVector<> getLowerBounds(gsVector<> upperBounds);

  void getBounds(patchSide ps, gsVector<> &bounds, index_t fixInfo);

  IpOptSparseMatrix getJacobian() { return m_constraintMatrix; }

  void freeBoundary(index_t p, index_t bnd, index_t k);

  index_t findBoundaryIndex(index_t p, index_t bnd);


  friend class liaoOptProblem;

protected:
	gsMultiPatch<>* mp;

  size_t n_controlpoints;

  IpOptSparseMatrix m_constraintMatrix;

  gsVector<index_t> bndFixInfo;

public:
  real_t aBigNumber;
  index_t n_constraints;

};




#endif //INTERFACECONSTRAINT_H
