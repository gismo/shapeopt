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
    gsOptParamMethod(gsMultiPatch<>* mpin, bool use_dJC);

    // Constructs from mappers, one for each coordinate, should be finalized with
    // design variables for shape optimization tagged.
    gsOptParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers);

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
    // FIXIT use FD as default
    gsVector<> gradObj() const;

    // Evaluation of hessian of objective. Should be overloaded in inherited class
    // FIXIT use FD as default
    gsMatrix<> hessObj() const {};

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

    // Method to print optimization parameters
    void print();

public:
    bool use_detJacConstraint = false;
    mutable gsDetJacConstraint m_dJC;

    // FIXIT: should we add control over n.o. quadrature points, e.g. quA and quB?
};

#endif //GSOPTPARAMMETHOD_H
