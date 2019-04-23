#ifndef springMethod_H
#define springMethod_H
using namespace gismo;

#include "detJacConstraint.h"
#include "paraWithMapper.h"
#include "paraOptProblem.h"

class springMethod: public paraWithMapper{
public:
  springMethod(gsMultiPatch<>* mpin);
  springMethod(paraOptProblem* problem);

  // Methods for solving
  gsVector<> solveSystems(gsVector<> deltaCps);
  gsVector<> solve(gsVector<> deltaCps);
  void solve();
  void solveAndUpdate(gsVector<> deltaCps);

  index_t numDesignVars(){ return ndesign; }
  gsVector<> refControlPoints(){ return refCps; }

  void setupSystem(gsDofMapper &mapper, gsSparseMatrix<> &A, gsVector<> &b, index_t coord);
  bool is_boundary(index_t i, index_t p) const;

  // temporary method ...
  void reset(){};

protected:
  detJacConstraint dJC;
  // gsMultiPatch<>* mp is already in the base class paraWithMapper
  static const index_t d = 2; // Dimension 2

  gsVector<> refCps;

  std::vector<gsSparseMatrix<>> A; // Systems to find xcoordinates
  std::vector<gsVector<>> b; // Rhs
  index_t ndesign;
  std::vector< Eigen::FullPivLU<Eigen::Matrix<real_t,Dynamic,Dynamic>> > solvers;
};


#endif //springMethod_H
