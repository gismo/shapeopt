/** @file gsOptInit.h

@brief Creates initial guess.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/

#ifndef GSOPTINIT_H
#define GSOPTINIT_H
using namespace gismo;

#include "gsParamMethod.h"
#include "gsIpOptSparseMatrix.h"
#include <gsIpopt/gsOptProblem.h>

class gsOptInit: public gsParamMethod, public gsOptProblem<real_t>{
public:
    gsOptInit(memory::shared_ptr<gsMultiPatch<>> mpin, gsVector<> &tagged_goal);

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    void setupOptParameters();

    // Update controlpoints. Basicly calls the solve method from gsOptProblem<>
    // Returns false if the update failed
    bool update(){ solve(); return m_status; };

	gsVector<> evalCon() const;
	
	gsIpOptSparseMatrix jacobCon() const;

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

    // Method that is called between each optimization iteration.. Can be used to log
    //      or check stuff. If it returns false the optimization will be interrupted.
    bool intermediateCallback();

	void computeConstraintMatrix();

public:
	index_t n_boundaries;
	index_t n_tagged_cps;

	gsVector<index_t> m_global_to_tagged;

	gsVector<> m_tagged_goal;

	gsMatrix<> m_A; // Constraints
	gsIpOptSparseMatrix::Ptr m_A_ipopt; // Constraints in ipopt format.
};

#endif //GSOPTINIT_H
