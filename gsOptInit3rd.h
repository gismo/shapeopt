/** @file gsOptInit3rd.h

@brief Creates initial guess.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/

#ifndef GSOPTINIT3rd_H
#define GSOPTINIT3rd_H
using namespace gismo;

#include "gsParamMethod.h"
#include "gsIpOptSparseMatrix.h"
#include <gsIpopt/gsOptProblem.h>

class gsOptInit3rd: public gsParamMethod, public gsOptProblem<real_t>{
public:
    gsOptInit3rd(memory::shared_ptr<gsMultiPatch<>> mpin);

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    void setupOptParameters();

    // Update controlpoints. Basicly calls the solve method from gsOptProblem<>
    // Returns false if the update failed
    bool update(){ solve(); return m_status; };

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

	void getCorners();

	void updateCorners(gsMultiPatch<>::Ptr mp);

	bool isCorner(index_t ii) const;

	void markBoundaryAndOppositeToDof(patchSide &ps, gsMatrix< bool > &ps_flag);

	std::vector< patchSide > getUnmarkedNeighbours(patchSide &ps, gsMatrix< bool > &ps_flag);

	void markBoundaryToDof(boxSide &bs, index_t p, index_t dof, gsMatrix< bool > &ps_flag);

	void sort(gsMatrix< unsigned > &v);

	index_t countIntersection(gsMatrix< unsigned > &v, gsMatrix< unsigned > &w);


public:
	index_t n_boundaries;
	index_t n_boxes;
	index_t n_corners;
	index_t m_d;

	gsVector<index_t> m_count;

	gsVector<index_t> m_direction_map;

	gsMatrix<index_t> m_boundary_to_dof;

	gsVector<> m_goal;
};

#endif //GSOPTINIT3rd_H
