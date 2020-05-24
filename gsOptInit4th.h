/** @file gsOptInit4th.h

@brief Creates initial guess.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/

#ifndef GSOPTINIT4th_H
#define GSOPTINIT4th_H
using namespace gismo;

#include "gsParamMethod.h"
#include "gsIpOptSparseMatrix.h"
#include <gsIpopt/gsOptProblem.h>

class gsOptInit4th: public gsParamMethod, public gsOptProblem<real_t>{
public:
    gsOptInit4th(gsMultiPatch<>::Ptr mpin, gsMultiPatch<>::Ptr mp_goal);

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    void setupOptParameters();

    // Update controlpoints. Basicly calls the solve method from gsOptProblem<>
    // Returns false if the update failed
    bool update(){ solve(); return m_status; };

    real_t evalObj ( ) const;
    gsVector<> gradObj ( ) const;

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

    void print();

    // Method that is called between each optimization iteration.. Can be used to log
    //      or check stuff. If it returns false the optimization will be interrupted.
    bool intermediateCallback();

    index_t getDesignIndex(const index_t i, const index_t j) const;
    index_t getDesignIndex(const index_t i) const;

    void updateFromDesignVars( const gsVector<> & u) const;
    gsVector<> getDesignVars() const;

    void updateMP(); // update m_mp from m_A and m_b

    gsMatrix<> A(){ return m_A; };
    gsMatrix<> b(){ return m_b; };

public:
	index_t m_d;
    index_t n_patches;
    index_t n_boundaries;
    index_t n_design;

    mutable gsMatrix<> m_A;
    mutable gsVector<> m_b;

    gsMultiPatch<>::Ptr m_mp_goal;
};

#endif //GSOPTINIT4th_H
