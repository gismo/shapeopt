#ifndef linearizedOptProblem_H
#define linearizedOptProblem_H
using namespace gismo;

#include "paraOptProblem.h"

class linearizedOptProblem{
public:
  linearizedOptProblem(paraOptProblem* problem, gsMatrix<> constMat,gsVector<> constVals):
      m_problem(problem), m_cMatrix(constMat), m_cValues(constVals){};
  linearizedOptProblem(paraOptProblem* problem);
  gsVector<> solve(gsVector<> deltaCps);
  void solveAndUpdate(gsVector<> deltaCps);
  real_t evalObj(real_t obj, gsVector<> grad, gsMatrix<> hess, gsVector<> x, gsVector<> x0);

  index_t numDesignVars(){ return ndesign; }
  gsVector<> refControlPoints(){ return refCps; }

protected:
  paraOptProblem* m_problem;
  gsMatrix<> m_cMatrix;
  gsVector<> m_cValues;

  gsVector<> refCps;
  gsVector<> grad;
  gsMatrix<> hess;
  real_t obj;

  gsMatrix<> M;
  gsMatrix<> KKTsystem;
  // gsSparseMatrix<> KKTsystem;
  gsVector<> rhs;
  index_t nconst;
  index_t ndesign;
  index_t nz;
  gsVector<> diffInBounds;
  Eigen::FullPivLU<Eigen::Matrix<real_t,Dynamic,Dynamic>> solver;
	// gsSparseSolver<>::LU solver;
};


#endif //linearizedOptProblem_H
