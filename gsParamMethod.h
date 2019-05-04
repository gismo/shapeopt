/** @file gsIpOptSparseMatrix.h

@brief  A Implements a base class for parametrization method for multipatch geometry
        It uses a list of gsDofMappers (one for each coordinate) which to denote wich
        DoFs that are fixed and which are not fixed.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSPARAMMETHOD_H
#define GSPARAMMETHOD_H
using namespace gismo;

class gsParamMethod{
public:

  // Constructs from multipatch, should setup mappers and assume interfaces glued
  // and boundary fixed; FIXIT: implement this functionality, FIXIT: tag bnd?
  gsParamMethod(gsMultiPatch<>* mpin);

  // Constructs from list of mappers
  gsParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers);

  // FIXIT: Think about which one of these functions is the best to have virtual...
  //        i.e. one could be implemented generic by the other one...
  // Method to update inner controlpoints. Should be implemented in inherited classes
  virtual void update() = 0;

  // Method to update inner controlpoints, given tagged Dofs (x).
  // Should be implemented in inherited classes
  virtual void update(gsVector<> x) = 0;

  // Method to get the jacobian of the update (free Dofs) with respect to x (the tagged Dofs)
  // Should be implemented in inherited classes
  virtual gsMatrix<> jacobUpdate(gsVector<> x) = 0;

  // Get a vector of the tagged control points.
  gsVector<> getTagged();

  // Update the tagged control points.
  void updateTagged(gsVector<> x) const;

  // Returns vector with free variables (with coordinates stacked, first x then y the possibly z)
  gsVector<> getFree() const;

  // Update the gsMultiPatch* mp, based on a vector of free variables
  void updateFree(gsVector<> des) const;

  // Update the gsMultiPatch* mp, based on a vector of free variables, and vector of tagged cps
  void updateFreeAndTagged(gsVector<> des, gsVector<> x) const;

  // Get the vector of all control points (eliminated and free)
  gsVector<> getControlPoints() const;

  // Get the vector of all control points (eliminated and free) from design variables
  gsVector<> getControlPoints(gsVector<> des) const;

  // Updates m_mp by vector of all controlpoints (e.g. obtained from getControlPoints() )
  // OBS: also updates eliminated DoFs...
  void updateControlPoints(gsVector<> cps);

  // Mappers for antenna problem, FIXIT: should be moved to gsOptAntenna at some point.
  void setupMapper();


public:
  mutable gsMultiPatch<>* m_mp;

  std::vector< gsDofMapper > m_mappers; // Mapper for each coordinate

  gsVector<> m_shift_free; // Contains the shifts used when mapping free DoFs
  gsVector<> m_shift_all; // Contains the shifts used when mapping all DoFs

  index_t n_free; // sum over coordinates
  index_t n_cps; // sum over coordinates
  index_t n_tagged; // sum over coordinates
  index_t n_controlpoints;
  index_t fixedPatch = 3; // FIXIT: Should be moved to gsOptAntenna.

};

# endif //GSPARAMMETHOD_H
