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
#include "gsShapeOptLog.h"
#include <gsIpopt/gsOptProblem.h>

class gsOptParamMethod: public gsParamMethod, public gsOptProblem<real_t>{
public:
    // Constructs from multipatch, by eliminating boundary, gluing interfaces and tagging boundary.
    // Uses gsDetJacConstraint or no constraints depending on \a use_dJC.
    // The gsDetJacConstraint exploinds tensor product structure depending on \a useTensorStructureforDJC.
    gsOptParamMethod(memory::shared_ptr<gsMultiPatch<>> mpin, bool use_dJC = true, bool useTensorStructureforDJC = false);

    // Constructs from mappers, one for each coordinate, should be finalized with
    // design variables for shape optimization tagged.
    // Uses gsDetJacConstraint or no constraints depending on \a use_dJC.
    // The gsDetJacConstraint exploinds tensor product structure depending on \a useTensorStructureforDJC.
    gsOptParamMethod(memory::shared_ptr<gsMultiPatch<>> mpin, std::vector< gsDofMapper > mappers, bool use_dJC = true, bool useTensorStructureforDJC = false);

    // Constructs from multipatch, by eliminating boundary, gluing interfaces and tagging boundary.
    // Uses the constraint given as input \a constraint
    gsOptParamMethod(memory::shared_ptr<gsMultiPatch<>> mpin, memory::shared_ptr<gsConstraint> constraint);

    // Constructs from mappers, one for each coordinate, should be finalized with
    // design variables for shape optimization tagged.
    gsOptParamMethod(memory::shared_ptr<gsMultiPatch<>> mpin, std::vector< gsDofMapper > mappers, memory::shared_ptr<gsConstraint> constraint);

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    void setupOptParameters();

    // Update controlpoints. Basicly calls the solve method from gsOptProblem<>
    // Returns false if the update failed
    bool update();

    // Update controlpoints from tagged. Left empty for now..
    // Returns false if the update failed
    virtual bool update(gsVector<> x){ GISMO_ERROR("Not implemented!");  };

    // Derivatives of update with respect to tagged Dofs. Left empty for now.
    virtual gsMatrix<> jacobUpdate(gsVector<> x){ GISMO_ERROR("Not implemented!");  };

    // Evaluation of objective. Should be overloaded in inherited class
    virtual real_t evalObj () const = 0;

    // Evaluation of gradient of objective.
    // Uses FD as default
    // Should be overloaded in inherited class, if exact gradient is achievable
    virtual gsVector<> gradObj() const;

    // Evaluation of hessian of objective. Calls hessAll;
    // Returns hessian wrt. free variables (Hcc), and hessian wrt to tagged and free in
    // \a Hcx.
    // FIXIT use FD as default
    virtual gsMatrix<> hessObj(gsMatrix<> &Hcx) const;

    // Evaluation of hessian of objective wrt all controlpoints
    // SHould be overloaded in inherited class.
    virtual gsMatrix<> hessAll(gsDofMapper &space_mapper) const
    {
        GISMO_NO_IMPLEMENTATION;
    };

    // Evaluation of jacobian of constraints. Calls mapMatrix.
    gsIpOptSparseMatrix jacobCon() const;

    //  Method to evaluate Lagrangian, uses lagrange multipliers from m_lambda (see gsOptProblem.h)
    real_t evalLagrangian () const;

    //  Method to evaluate gradient of Lagrangian, uses lagrange multipliers from m_lambda (see gsOptProblem.h)
    gsVector<> gradLagrangian () const;

    //  Method to evaluate hessian Lagrangian, uses lagrange multipliers from m_lambda (see gsOptProblem.h)
    gsMatrix<> hessLagrangian(gsMatrix<> &Hcx) const;

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
    void refineBasedOnDetJ(index_t strategy, memory::shared_ptr<gsDetJacConstraint> dJC);

    void setLog(memory::shared_ptr<gsShapeOptLog> inLog, memory::shared_ptr<gsDetJacConstraint> d){ m_d = d; m_log = inLog; use_log = true; };

    // Method that is called between each optimization iteration.. Can be used to log
    //      or check stuff. If it returns false the optimization will be interrupted.
    bool intermediateCallback();

    void setIntegrationBasis(gsMultiBasis<> &mesh) { m_integrationBasis = memory::make_shared(&mesh); };

public:
    bool use_detJacConstraint = false;
    memory::shared_ptr<gsConstraint> m_dJC;

    memory::shared_ptr<gsShapeOptLog> m_log;
    bool use_log = false;
    memory::shared_ptr<gsDetJacConstraint> m_d; //tmp

    index_t m_iter = 0; // Iteration count

    // Control over no quadrature points. Set by setQuad(quA,quB) method.
    real_t m_quA = 4;
    real_t m_quB = 4;

    // pointer to basis on which to perform integration
    // It is derived from *m_mp as default, but can be set to some other basis
    // This could be a locally refined mesh
    memory::shared_ptr<gsMultiBasis<>> m_integrationBasis;
};

#endif //GSOPTPARAMMETHOD_H
