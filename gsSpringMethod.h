/** @file gsSpringMethod.h

@brief  Implements the spring parametrization method. Inherits from gsAffineParamMethod.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef gsSpringMethod_H
#define gsSpringMethod_H
using namespace gismo;

#include "gsAffineParamMethod.h"

class gsSpringMethod: public gsAffineParamMethod{
public:
  // Constructs using gsMultiPatch, by eliminating boundaries and gluing interfaces
  gsSpringMethod(gsMultiPatch<>* mpin);

  // Constructs using mappers. There should be one mapper for each coordinate.
  // And they should already be finalized and the firstIndex of the second two
  // should be shifted with the freeSize of the previous ones
  gsSpringMethod(gsMultiPatch<>* mpin,std::vector< gsDofMapper > mappers);

  // Returns the free Dofs from tagged cps
  gsVector<> getUpdate(gsVector<> x);

  // setup up system from one mapper
  void setupSystem(gsDofMapper &mapper, gsSparseMatrix<> &A, gsVector<> &b, index_t coord);

  // Setup up all d systems using setupSystem for each of them
  void setupSolvers();

  // Method to find it local Dof is on boundary
  // FIXIT: very expensive, find better method
  bool is_boundary(index_t i, index_t p) const;

  // temporary method ...
  void reset(){};

protected:
  // FIXIT: make dimension independent
  static const index_t d = 2; // Dimension 2

  // FIXIT: do I need to save A and b or is solvers enough?
  std::vector<gsSparseMatrix<>> m_As; // Systems to find coordinates
  std::vector<gsVector<>> m_bs; // Rhs
  std::vector< Eigen::FullPivLU<Eigen::Matrix<real_t,Dynamic,Dynamic>> > m_solvers;
};


#endif //gsSpringMethod_H
