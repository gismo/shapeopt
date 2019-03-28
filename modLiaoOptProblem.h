#ifndef modLiaoOptProblem_H
#define modLiaoOptProblem_H
using namespace gismo;

#include "detJacConstraint.h"
#include "interfaceConstraint.h"
#include "paraOptProblem.h"


class modLiaoOptProblem: public paraOptProblem {
public:
  modLiaoOptProblem(gsMultiPatch<>* mpin):paraOptProblem(mpin){};
  real_t evaluateOnPatch(index_t i) const;
  void evaluateDerivOnPatch(index_t i, gsVector<> &xVec, gsVector<> &yVec) const;
  void evaluate2ndDerivOnPatch(index_t i, gsMatrix<> &xxMat, gsMatrix<> &xyMat, gsMatrix<> &yyMat) const;

  gsVector<> gradObjFD(gsMatrix<> u) const;
  void setNoQuadraturePoints(real_t a, int b){
    quA = a; quB = b;
  }
  real_t quA = 1;
  int quB = 1;
};

#endif
