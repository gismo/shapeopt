#ifndef paraOptProblemBase_H
#define paraOptProblemBase_H
using namespace gismo;

#include <gsIpopt/gsOptProblem.h>

class paraOptProblemBase: public gsOptProblem<real_t> {
public:
  virtual real_t evalObj () const = 0;
  virtual gsVector<> gradientObj() const = 0;
  virtual gsMatrix<> hessianObj() const = 0;

  virtual gsVector<> getDesignVariables() const = 0;
  virtual void updateDesignVariables(gsVector<> des) = 0;

  int numConJacNonZero(){ return m_numConJacNonZero; }
  gsVector<> desLowerBounds() const { return m_desLowerBounds; }
  gsVector<> desUpperBounds() const { return m_desUpperBounds; }


};

# endif //paraOptProblemBase_H
