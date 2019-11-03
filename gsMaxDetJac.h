#ifndef GSMAXDETJAC_H
#define GSMAXDETJAC_H
using namespace gismo;

#include <gsIpopt/gsOptProblem.h>
#include "gsParamMethod.h"
#include "gsDetJacConstraint.h"
#include "gsIpOptSparseMatrix.h"

class gsMaxDetJac:  public gsParamMethod, public gsOptProblem<real_t>
{
public:

    // Constructs from multipatch, by eliminating boundary, gluing interfaces and tagging bnd
    gsMaxDetJac(gsMultiPatch<>* mpin, bool useTPSolver = false);

    // Constructs from mappers, one for each coordinate, should be finalized with
    // design variables for shape optimization tagged.
    gsMaxDetJac(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers, bool useTPSolver = false);

    // Method to set the optimization parameters such as design bounds, constraint bounds etc   .
    void setupOptParameters();

    // Update controlpoints. Basicly calls the solve method from gsOptProblem<>
    // Returns false if the update failed
    bool update();

    // Update controlpoints from tagged. Left empty for now..
    // Returns false if the update failed
    bool update(gsVector<> x){ GISMO_ERROR("Not implemented!");  };

    // Derivatives of update with respect to tagged Dofs. Left empty for now.
    gsMatrix<> jacobUpdate(gsVector<> x){ GISMO_ERROR("Not implemented!");  };

    real_t evalObj ( const gsAsConstVector<real_t> & u) const;
    void gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const;
    void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;
    void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const ;

    // Method overloaded from gsOptProblem<>
    // Uses gsDetJacConstraint depending on the flag use_detJacConstraint
    void computeJacStructure();

    // Print the no design variable, constraints etc.
    void print();


public:
    mutable gsDetJacConstraint m_dJC;

};

# endif //GSMAXDETJAC_H
