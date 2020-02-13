#include <gismo.h>
#include "gsMaxDetJac.h"
using namespace gismo;

gsMaxDetJac::gsMaxDetJac(memory::shared_ptr<gsMultiPatch<>> mpin, bool useTPSolver):
gsParamMethod(mpin), m_dJC(mpin, useTPSolver)
{
    setupOptParameters();
};

gsMaxDetJac::gsMaxDetJac(memory::shared_ptr<gsMultiPatch<>> mpin, std::vector< gsDofMapper > mappers, bool useTPSolver):
gsParamMethod(mpin, mappers), m_dJC(mpin, useTPSolver)
{
    setupOptParameters();
};

bool gsMaxDetJac::update()
{
    solve();

    return true; // FIXIT: return status IpOpt::Solve_Succeeded
};

void gsMaxDetJac::setupOptParameters()
{
    // The desing variables is the free variables and a slack variable
    m_curDesign.setZero(n_free + 1,1);
    m_curDesign.block(0,0,n_free,1) = getFree();
    m_curDesign(n_free,0) = m_dJC.evalCon().minCoeff();

    m_numDesignVars = n_free + 1;

    // Essentially no design bounds
    m_desLowerBounds.setConstant(n_free + 1, -1e9);
    m_desUpperBounds.setConstant(n_free + 1,  1e9);

    // Constraint bounds is given by gsDetJacConstraint
    // Call once to setup solver
    m_dJC.evalCon();

    m_numConstraints = m_dJC.numConstraints();
    m_conLowerBounds.setZero(m_numConstraints);
    m_conUpperBounds.setConstant(m_numConstraints,1e9);

    // compute jac structure
    computeJacStructure();
}

real_t gsMaxDetJac::evalObj ( const gsAsConstVector<real_t> & u) const {
    return -u[m_numDesignVars - 1]; // Return slack variable
}

void gsMaxDetJac::gradObj_into(const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const{
    // Set almost almost the most
    for (index_t i = 0; i < m_numDesignVars - 1; i++){
        result[i] = 0;
    }
    result[m_numDesignVars - 1] = -1;
}

void gsMaxDetJac::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
    // gsInfo << "...evalCon_into" << "\n" << std::flush;

    gsVector<> des = u.segment(0,m_numDesignVars-1);
    updateFree(des);

    // we consider constraints d_i > slack for all i
    // OBS: assumes positive determinant
    result.segment(0,m_dJC.n_constraints) = m_dJC.evalCon();

    for (index_t i = 0; i < m_dJC.n_constraints; i++){
        result[i] -= u[m_numDesignVars - 1];
    }

    // gsInfo << "Max of constraint " << result.segment(0,dJC.n_constraints).maxCoeff() << "\n";
}

void gsMaxDetJac::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
    // gsInfo << "...jacobCon_into" << "\n" << std::flush;
    updateFree(u.segment(0,m_numDesignVars-1));

    // Compute the first part, namely the gradient of the det Jac constraints
    gsIpOptSparseMatrix J = mapMatrix(m_dJC.space_mapper(),m_dJC.getJacobian());

    // Compute final part of gradient
    gsMatrix<> ones;
    ones.setConstant(m_numConstraints,1,-1);
    gsIpOptSparseMatrix J2(ones,0);

    // Concatenate J and J2
    J.concatenate(J2,"row");

    // gsInfo << "We expect nnz " << m_numConJacNonZero << "\n";
    // gsInfo << "We have in J1 " << J1.nnz() << "\n";
    result = J.values();
}

void gsMaxDetJac::computeJacStructure()
{
    // Use the sparsity provided by gsDetJacConstraint m_dJC
    gsIpOptSparseMatrix J = mapMatrix(m_dJC.space_mapper(),m_dJC.getJacobian());

    // Compute final part of gradient
    gsMatrix<> ones;
    ones.setConstant(m_numConstraints,1,-1);
    gsIpOptSparseMatrix J2(ones,0);

    // Concatenate J and J2
    J.concatenate(J2,"row");

    m_numConJacNonZero = J.nnz();
    m_conJacRows = J.rows();
    m_conJacCols = J.cols();

};

void gsMaxDetJac::print(){
  gsInfo << "m_numDesignVars  = " <<  m_numDesignVars << "\n";
  gsInfo << "m_numConstraints  = " <<  m_numConstraints << "\n";

  gsInfo << "m_numConJacNonZero  = " <<  m_numConJacNonZero << "\n";
  gsInfo << "m_conJacRows.size  = " <<  m_conJacRows.size() << "\n";

  // for(std::vector<index_t>::iterator it = m_conJacRows.begin(); it != m_conJacCols.end(); ++it){
  //     gsInfo << " " << *it;
  // }

  gsInfo << "\nm_conJacCols.size  = " <<  m_conJacCols.size() << "\n";

  gsInfo << "m_conUpperBounds.size() = " << m_conUpperBounds.size() << "\n";
  gsInfo << "m_conLowerBounds.size() = " << m_conLowerBounds.size() << "\n";
  gsInfo << "m_desUpperBounds.size() = " << m_desUpperBounds.size() << "\n";
  gsInfo << "m_desLowerBounds.size() = " << m_desLowerBounds.size() << "\n";

  gsInfo << "m_curDesign.size() = " << m_curDesign.size() << "\n";

  gsMatrix<> disp(m_desUpperBounds.size(),2);
  disp << m_desLowerBounds,m_desUpperBounds;
  // gsInfo << ".. design upper and lower bounds\n";
  // gsInfo << disp << "\n";
  //
  gsMatrix<> disp2(m_conUpperBounds.size(),2);
  disp2 << m_conLowerBounds,m_conUpperBounds;
  // gsInfo << ".. constraint upper and lower bounds\n";
  // gsInfo << disp2 << "\n";

}
