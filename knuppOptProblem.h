#ifndef knuppOptProblem_H
#define knuppOptProblem_H
using namespace gismo;

#include "detJacConstraint.h"
#include "interfaceConstraint.h"
#include "paraOptProblem.h"


class knuppOptProblem: public paraOptProblem {
public:
  knuppOptProblem(gsMultiPatch<>* mpin):paraOptProblem(mpin){};
  real_t evaluateOnPatch(index_t i) const;
  void evaluateDerivOnPatch(index_t i, gsVector<> &xVec, gsVector<> &yVec) const;
  void evaluate2ndDerivOnPatch(index_t i, gsMatrix<> &xxMat, gsMatrix<> &xyMat, gsMatrix<> &yyMat) const;

private:
  real_t k_A = 0.5; //Weight on Area part of Knupps functional
  real_t k_O = 0.5; //Weight on Orthogonality part of Knupps functional
};

#endif //knuppOptProblem_H
