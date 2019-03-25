#ifndef paraOptProblem_H
#define paraOptProblem_H
using namespace gismo;

#include "detJacConstraint.h"
#include "interfaceConstraint.h"
#include "paraOptProblemBase.h"


class paraOptProblem: public paraOptProblemBase {
public:
  paraOptProblem(gsMultiPatch<>* mpin);
  gsVector<> getDesignVariables() const;
  void updateDesignVariables(gsVector<> des);
  real_t evalObj() const;
  void gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const;
  real_t evalObj ( const gsAsConstVector<real_t> & u) const;
  void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;
  void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;
  void print();
  int numConJacNonZero(){ return m_numConJacNonZero; }
  void writeToFile(gsVector<> vec, std::string name) const;
  void writeToFile(gsMatrix<> mat, std::string name) const;
  void loadFromFile(std::string name);
  gsMatrix<> hessianObj() const;
  gsVector<> gradientObj() const;

  //Virtual methods to implement in derived classes
  virtual real_t evaluateOnPatch(index_t i) const = 0;
  virtual void evaluateDerivOnPatch(index_t i, gsVector<> &xVec, gsVector<> &yVec) const = 0;
  virtual void evaluate2ndDerivOnPatch(index_t i, gsMatrix<> &xxMat, gsMatrix<> &xyMat, gsMatrix<> &yyMat) const = 0;

  gsVector<> desLowerBounds() const { return m_desLowerBounds; }
  gsVector<> desUpperBounds() const { return m_desUpperBounds; }


public:
  mutable gsMultiPatch<>* mp;
  mutable detJacConstraint dJC;
  mutable interfaceConstraint iC;
  gsMatrix<> interfaceConstraintMatrix;
  real_t m_eps = 0.0001;
  mutable index_t counter1 = -2;
  mutable index_t counter2 = -2;
  index_t fixedPatch = 3;


};

#endif //paraOptProblem_H
