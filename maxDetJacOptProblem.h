#ifndef maxDetJacOptProblem_H
#define maxDetJacOptProblem_H
using namespace gismo;

#include <gsIpopt/gsOptProblem.h>
#include "detJacConstraint.h"
#include "interfaceConstraint.h"
#include "paraOptProblem.h"

class maxDetJacOptProblem: public paraOptProblem {
public:
  maxDetJacOptProblem(gsMultiPatch<>* mpin);

  real_t evalObj ( const gsAsConstVector<real_t> & u) const;
  void gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const;
  void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;
  void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

  // Methods to reset the fixed degrees of freedoms
  void setDesignBounds();
  void reset();                 //calls setDesignBounds() and set m_curDesign to zero;

  void writeToFile(gsVector<> vec, std::string name) const;
  void writeToFile(gsMatrix<> mat, std::string name) const;

  // Has to be implemented unfortunately FIXIT TODO
  real_t evaluateOnPatch(index_t i) const{return 0;};
  void evaluateDerivOnPatch(index_t i, gsVector<> &xVec, gsVector<> &yVec) const{};
  void evaluate2ndDerivOnPatch(index_t i, gsMatrix<> &xxMat, gsMatrix<> &xyMat, gsMatrix<> &yyMat) const{};


public:
  mutable gsMultiPatch<>* mp;
  mutable detJacConstraint dJC;
  mutable interfaceConstraint iC;
  gsMatrix<> interfaceConstraintMatrix;
  index_t fixedPatch = 3;
};

# endif //paraOptProblemBase_H
