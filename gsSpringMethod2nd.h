/** @file gsSpringMethod2nd.h

@brief  Implements the spring parametrization method. Inherits from gsAffineParamMethod. It is altered such that bnd cps should be an average of its boundary neighbors.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/

#ifndef gsSpringMethod2nd_H
#define gsSpringMethod2nd_H
using namespace gismo;

#include "gsSpringMethod.h"

class gsSpringMethod2nd: public gsSpringMethod{
public:
  // Copy constuctors from Spring method
  using gsSpringMethod::gsSpringMethod;

  // setup up system from one mapper
  void setupSystem(gsDofMapper &mapper, gsSparseMatrix<> &A, gsVector<> &b, index_t coord);

  void getDirections(std::vector< index_t > &result, index_t i, index_t p);

  bool is_double_boundary(index_t i, index_t p);

};


#endif //gsSpringMethod2nd_H
