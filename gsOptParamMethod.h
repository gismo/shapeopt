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
#include "gsIpOptSparseMatrix.h"
#include <gsIpopt/gsOptProblem.h>

class gsOptParamMethod: public gsParamMethod, public gsOptProblem<real_t>{
public:
    // Constructs from multipatch, by eliminating boundary, gluing interfaces and tagging bnd
    gsOptParamMethod(gsMultiPatch<>* mpin, bool use_dJC = true, bool useTensorStructureforDJC = false);

    // Constructs from mappers, one for each coordinate, should be finalized with
    // design variables for shape optimization tagged.
    gsOptParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers, bool use_dJC, bool useTensorStructureforDJC = false);

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

    // Evaluation of jacobian of constraints. Calls mapMatrix.
    gsIpOptSparseMatrix jacobCon() const;

    //  Method to evaluate Lagrangian, uses lagrange multipliers from m_lambda (see gsOptProblem.h)
    real_t evalLagrangian () const;

    //  Method to evaluate gradient of Lagrangian, uses lagrange multipliers from m_lambda (see gsOptProblem.h)
    gsVector<> gradLagrangian () const;

    //  Method to evaluate hessian Lagrangian, uses lagrange multipliers from m_lambda (see gsOptProblem.h)
    gsMatrix<> hessLagrangian(gsMatrix<> &hessObjTagged) const;

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

    // set quadrature parameters
    void setQuad(real_t quA, index_t quB){ m_quA = quA, m_quB = quB; };

    // Refine elements OBS: requires mp to have gsHTensorBasis
    // Strategy can be:
    //      0 - refine where detJ is negative,
    //      1 - refine where detJ constraints are active
    // Remember that you have to update mappers afterwards
    void refineBasedOnDetJ(index_t strategy);

public:
    bool use_detJacConstraint = false;
    mutable gsDetJacConstraint m_dJC;

    // Control over no quadrature points. Set by setQuad(quA,quB) method.
    real_t m_quA = 2;
    real_t m_quB = 2;
};

#endif //GSOPTPARAMMETHOD_H
