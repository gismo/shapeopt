#ifndef harmonicOptProblem_H
#define harmonicOptProblem_H
using namespace gismo;

#include "detJacConstraint.h"
#include "interfaceConstraint.h"
#include <gsIpopt/gsOptProblem.h>

class harmonicOptProblem: public gsOptProblem<real_t> {
public:
  harmonicOptProblem(gsMultiPatch<>* mpin);
  real_t evalObj () const;
  // virtual gsVector<> gradientObj() const = 0;
  // virtual gsMatrix<> hessianObj() const = 0;

  gsVector<> getDesignVariables() const;
  void updateDesignVariables(gsVector<> des);

  int numConJacNonZero(){ return m_numConJacNonZero; }
  gsVector<> desLowerBounds() const { return m_desLowerBounds; }
  gsVector<> desUpperBounds() const { return m_desUpperBounds; }

  real_t evalObj ( const gsAsConstVector<real_t> & u) const;
  void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;
  void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;
  void print();

  void setDesignBounds();
  void reset();                 //calls setDesignBounds() and set m_curDesign to zero;

public:
  mutable gsMultiPatch<>* mp;
  mutable detJacConstraint dJC;
  mutable interfaceConstraint iC;
  gsMatrix<> interfaceConstraintMatrix;
  real_t m_eps = 0.0001;
  index_t fixedPatch = 3;

  real_t lambda_1 = 0.05;
  real_t lambda_2 = 0.05;



};

# endif //harmonicOptProblem_H
