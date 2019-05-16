/** @file gsOptParamMethod.h

@brief  Base class for parametrization methods based on nonlinear optimization problems,
        inherits from gsParamMethod (cps handling) and gsOptProblem<> (optimization with IpOpt)

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/

#ifndef GSOPTPARAMMETHOD_H
#define GSOPTPARAMMETHOD_H
using namespace gismo;

#include "gsParamMethod.h"
#include "gsDetJacConstraint.h"
#include <gsIpopt/gsOptProblem.h>

class gsOptParamMethod: public gsParamMethod, public gsOptProblem<real_t>{
public:
    // Constructs from multipatch, by eliminating boundary, gluing interfaces and tagging bnd
    // FIXIT: implement functionality
    gsOptParamMethod(gsMultiPatch<>* mpin, bool use_dJC);

    // Constructs from mappers, one for each coordinate, should be finalized with
    // design variables for shape optimization tagged.
    gsOptParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers, bool use_dJC);

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    void setupOptParameters();

    // Update controlpoints. Basicly calls the solve method from gsOptProblem<>
    void update();

    // Update controlpoints from tagged. Left empty for now..
    void update(gsVector<> x){ GISMO_ERROR("Not implemented!");  };

    // Derivatives of update with respect to tagged Dofs. Left empty for now.
    gsMatrix<> jacobUpdate(gsVector<> x){ GISMO_ERROR("Not implemented!");  };

    // Evaluation of objective. Should be overloaded in inherited class
    virtual real_t evalObj () const = 0;

    // Evaluation of gradient of objective.
    // Uses FD as default
    // Should be overloaded in inherited class, if exact gradient is achievable
    virtual gsVector<> gradObj() const;

    // Evaluation of hessian of objective. Should be overloaded in inherited class
    // FIXIT use FD as default
    virtual gsMatrix<> hessObj(gsMatrix<> &hessObjTagged) const {};

    // Method overloaded from gsOptProblem<>
    real_t evalObj ( const gsAsConstVector<real_t> & u) const;

    // Method overloaded from gsOptProblem<>
    void gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const;

    // Method overloaded from gsOptProblem<>
    // Uses gsDetJacConstraint depending on the flag use_detJacConstraint
    void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    // Method overloaded from gsOptProblem<>
    // Uses gsDetJacConstraint depending on the flag use_detJacConstraint
    void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    // Method overloaded from gsOptProblem<>
    // Uses gsDetJacConstraint depending on the flag use_detJacConstraint
    void computeJacStructure();

    // Maps a vector from mapper_in indexing, to m_mappers.
    // E.g. used to map gradients with respect to all cps to be respect to free
    // DoFs...
    // FIXIT: maybe find a more elegant way to handle all of this..
    //      Maybe always assemble wrt. specific mapper, and hold a permutation
    //      in a matrix, that we only need to multiply to map...
    gsMatrix<> mapMatrix(gsDofMapper mapper_in, gsMatrix<> mat) const;
    gsMatrix<> mapMatrixToTagged(gsDofMapper mapper_in, gsMatrix<> mat) const;

    // Maps a gsIpOptSparseMatrix from mapper_in indexing, to m_mappers.
    // E.g. used to map gradients with respect to all cps to be respect to free
    // DoFs...
    gsIpOptSparseMatrix mapMatrix(gsDofMapper mapper_in, gsIpOptSparseMatrix M) const;

    // Method to print optimization parameters
    void print();

public:
    bool use_detJacConstraint = false;
    mutable gsDetJacConstraint m_dJC;

    // FIXIT: should we add control over n.o. quadrature points, e.g. quA and quB?
};

#endif //GSOPTPARAMMETHOD_H
