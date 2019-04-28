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
  // and boundary fixed; FIXIT: implement this functionality
  gsParamMethod(gsMultiPatch<>* mpin);

  // Constructs from list of mappers
  gsParamMethod(std::vector< gsDofMapper > mappers);

  // Method to update inner controlpoints. Should be implemented in inherited classes
  virtual void update() = 0;

  // Method to get the jacobian of the update. Should be implemented in inherited classes
  virtual void jacobUpdate() = 0;

  // Returns vector with free variables (with coordinates stacked, first x then y the possibly z)
  gsVector<> getDesignVariables() const;

  // Update the gsMultiPatch* mp, based on a vector of free variables
  void updateDesignVariables(gsVector<> des) const;

  // Get the full vector on control points from patch
  gsVector<> getControlPoints() const;

  // Get the full vector on control points from design variables
  gsVector<> getControlPoints(gsVector<> des) const;

  // Mappers for antenna problem, FIXIT: should be moved to gsOptAntenna at some point.
  void setupMapper();

public:
  mutable gsMultiPatch<>* m_mp;

  std::vector< gsDofMapper > m_mappers; // Mapper for each coordinate

  index_t n_freeCps;
  index_t n_controlpoints;
  index_t fixedPatch = 3; // FIXIT: Should be moved to gsOptAntenna.

};

# endif //GSPARAMMETHOD_H
