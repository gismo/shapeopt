#ifndef paraWithMapper_H
#define paraWithMapper_H
using namespace gismo;

class paraWithMapper{
public:
  paraWithMapper(gsMultiPatch<>* mpin);
  gsVector<> getDesignVariables() const;
  void updateDesignVariables(gsVector<> des) const;

  gsVector<> getControlPoints() const; // Get the full vector on control points from patch
  gsVector<> getControlPoints(gsVector<> des) const; // Get the full vector on control points from design variables

  void setupMapper();

public:
  mutable gsMultiPatch<>* mp;

  std::vector< gsDofMapper > m_mappers; // Mapper for each coordinate

  index_t n_freeCps;
  index_t n_controlpoints;
  index_t fixedPatch = 3;

};

# endif //paraWithMapper_H
