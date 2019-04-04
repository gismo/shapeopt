#ifndef SHAPEOPTPROBLEM_H
#define SHAPEOPTPROBLEM_H
using namespace gismo;

#include "detJacConstraint.h"
#include "interfaceConstraint.h"
#include "modLiaoOptProblem.h"
#include "maxDetJacOptProblem.h"
#include "winslowOptProblem.h"
#include "linearizedOptProblem.h"
#include "IpOptSparseMatrix.h"
#include "stateEquationAntenna.h"
#include <gsIpopt/gsOptProblem.h>

class shapeOptProblem: public gsOptProblem<real_t>{
public:
  shapeOptProblem(gsMultiPatch<>* mpin,index_t numRefine, std::string output, bool plotDesign,
    bool plotMagnitude, bool plotSolution, bool saveCps);
  real_t evalObj() const ;
  gsVector<> gradientObj() const ;
  real_t evalObj( const gsAsConstVector<real_t> & u ) const;
  gsVector<> evalCon() const ;
  void gradObj_into ( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result ) const;
  void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;
  void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

  gsVector<> evaluateDerivativeTerm1(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;
  gsVector<> evaluateDerivativeTerm2(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;
  gsVector<> gradientObjWithoutAdjoint() const;
  gsVector<> getObjDerivativeDu(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const;

  void resetParametrizationToReference();

  void setCurrentDesign(gsVector<> x);
  gsVector<> getDesignVariables() const;
  void updateDesignVariables(gsVector<> u) const;
  gsMatrix<> extractPatchRows(index_t p, gsMatrix<> in) const;
  gsVector<> getUpdateToCps(gsVector<> u) const;
  gsMatrix<> derivativeOfDesignUpdate() const;
  IpOptSparseMatrix derivativeOfDetJac() const;
  void writeToFile(gsVector<> vec, std::string name) const;
  void plotDesignInParaview(std::string name);

  real_t volumeOfPatch(index_t p) const;
  gsMatrix<> derivVolumeOfPatch(index_t p) const;

  // Resets the reference parametrization
  void updateReferenceParametrization();
  void setDesignBounds();

  // Method to run full optimization strategy (including reparametrization)
  void runOptimization(index_t maxiter);
  bool intermediateCallback();

public:
  gsFunctionExpr<> delta;
  gsFunctionExpr<> ddeltadx;
  gsFunctionExpr<> ddeltady;

  mutable gsMultiPatch<>* mp;
  mutable detJacConstraint dJC;
  mutable interfaceConstraint iC;
  mutable modLiaoOptProblem pOP;
  mutable maxDetJacOptProblem mOP;
  mutable linearizedOptProblem linOP;
  mutable stateEquationAntenna SE;
  real_t m_eps = 0.0001;
  index_t antennaPatch = 3;
  gsVector<index_t> designIndiciesLocal;
  gsVector<index_t> designIndiciesGlobal;

  real_t d_g = 57.5*1e-9;       // Minimum allowed gap between the two antennas 57.5 [nm]
  real_t d_bw = 776.25*1e-9;    // Width of bounding box (Around x=0)
  real_t d_bh = 546.25*1e-9;    // Height of bounding box

  mutable index_t counter1 = 0;

  // Some options regarding output
  bool m_plotDesign;
  bool m_plotMagnitude;
  bool m_plotSolution;
  bool m_saveCps;
  std::string m_output;


};



# endif //SHAPEOPTPROBLEM_H
