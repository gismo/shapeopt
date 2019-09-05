#include <gismo.h>
#include "gsOptParamMethod.h"
using namespace gismo;

gsOptParamMethod::gsOptParamMethod(gsMultiPatch<>* mpin, bool use_dJC, bool useTensorStructureforDJC):
    gsParamMethod(mpin), m_dJC(mpin, useTensorStructureforDJC), use_detJacConstraint(use_dJC)
{
    setupOptParameters();
};

gsOptParamMethod::gsOptParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers, bool use_dJC, bool useTensorStructureforDJC):
        gsParamMethod(mpin, mappers), m_dJC(mpin), use_detJacConstraint(use_dJC)
{
    setupOptParameters();
};

void gsOptParamMethod::setupOptParameters()
{
    // The desing variables is the free variables
    m_curDesign = getFree();
    m_numDesignVars = n_free;

    // Essentially no design bounds
    m_desLowerBounds.setConstant(n_free, -1e9);
    m_desUpperBounds.setConstant(n_free, 1e9);

    // Constraint bounds is given by gsDetJacConstraint
    if (use_detJacConstraint){
        // Call once to setup solver
        m_dJC.evalCon();

        m_numConstraints = m_dJC.numConstraints();
        gsDebugVar(m_numConstraints);
        m_conLowerBounds = m_dJC.getLowerBounds();
        m_conUpperBounds = m_dJC.getUpperBounds();

        // Set lagrange multipliers to zeros initially
        m_lambda.setZero(m_numConstraints,1);

    } else {
        m_numConstraints = 0;
    }

    // compute jac structure
    computeJacStructure();
}

void gsOptParamMethod::update()
{
    // Update by solving the nonlinear optimization problem
    solve();

    gsInfo << "\nMax Lagrange multipliers from parametrization " << m_lambda.maxCoeff() << "\n";
    gsInfo << "Min Lagrange multipliers from parametrization " << m_lambda.minCoeff() << "\n\n";
};

real_t gsOptParamMethod::evalObj ( const gsAsConstVector<real_t> & u) const
{
    // gsInfo << "evalObj\n" << std::flush;
    updateFree(u);
    return evalObj();
};

void gsOptParamMethod::gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "gradObj_into\n" << std::flush;
    updateFree(u);
    result = gradObj();
};

void gsOptParamMethod::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "evalCon_into\n" << std::flush;
    updateFree(u);
    if (use_detJacConstraint){
        m_dJC.evalCon_into(result);
    } else {
        return;
    }
};

gsIpOptSparseMatrix gsOptParamMethod::jacobCon() const
{
    return mapMatrix(m_dJC.space_mapper(),m_dJC.getJacobian());
}

void gsOptParamMethod::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "jacobCon_into\n" << std::flush;
    updateFree(u);
    if (use_detJacConstraint){
        gsIpOptSparseMatrix J = jacobCon();
        result = J.values();
    } else {
        return;
    }

};

void gsOptParamMethod::computeJacStructure()
{
    // Use the sparsity provided by gsDetJacConstraint m_dJC
    if (use_detJacConstraint){
        gsIpOptSparseMatrix J = jacobCon();

        m_numConJacNonZero = J.nnz();
        m_conJacRows = J.rows();
        m_conJacCols = J.cols();
    } else {
        m_numConJacNonZero = 0;
    }

};

gsVector<> gsOptParamMethod::gradObj() const{
    gsVector<> u = getFree();
    const index_t n = u.rows();
    //GISMO_ASSERT((index_t)m_numDesignVars == n*m, "Wrong design.");

    gsVector<> result(n);

    gsMatrix<real_t> uu = u;//copy
    gsAsVector<real_t> tmp(uu.data(), n);
    gsAsConstVector<real_t> ctmp(uu.data(), n);
    index_t c = 0;

    // for all partial derivatives (column-wise)
    for ( index_t i = 0; i!=n; i++ )
    {
        // to do: add m_desLowerBounds m_desUpperBounds check
        tmp[i]  += real_t(0.00001);
        const real_t e1 = this->evalObj(ctmp);
        tmp[i]   = u[i] + real_t(0.00002);
        const real_t e3 = this->evalObj(ctmp);
        tmp[i]   = u[i] - real_t(0.00001);
        const real_t e2 = this->evalObj(ctmp);
        tmp[i]   = u[i] - real_t(0.00002);
        const real_t e4 = this->evalObj(ctmp);
        tmp[i]   = u[i];
        result[c++]= ( 8 * (e1 - e2) + e4 - e3 ) / real_t(0.00012);
    }

    return result;
}

void gsOptParamMethod::print()
{
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

// Strategy can be:
//      0 - refine where detJ is negative,
//      1 - refine where detJ constraints are active
void gsOptParamMethod::refineBasedOnDetJ(index_t strategy)
{
    std::vector<bool> elMarked;

    if (strategy == 0) // Mark support of basis function with negative coefficient
    {
        m_dJC.markElements(elMarked,-1);
    } else { // Mark support of basis function with active coefficient (tol close to lower bound)
        real_t tol = 0.0001;
        m_dJC.markElements(elMarked,tol);
    }

    refineElements(elMarked); // Refine m_mp

    recreateMappers();      // Recreate m_mappers, see gsParamMethod.h for further information
    m_dJC.setup();          // Reset the gsDetJacConstraint
    setupOptParameters();   // Reset optimization parameters

}

real_t gsOptParamMethod::evalLagrangian () const
{
    if (use_detJacConstraint) {
        return evalObj() + (m_lambda.transpose()*m_dJC.evalCon())(0); // element (0) is accesssed to convert to double
    } else { // If there are no constraints
        return evalObj();
    }
}

gsVector<> gsOptParamMethod::gradLagrangian () const
{
    if (use_detJacConstraint) {
        // FIXIT: Exploit sparsity
        gsMatrix<> J = jacobCon().asDense();
        return gradObj() + J.transpose()*m_lambda;
    } else { // If there are no constraints
        return gradObj();
    }
}

gsMatrix<> gsOptParamMethod::hessLagrangian(gsMatrix<> &hessObjTagged) const
{
    if (use_detJacConstraint) {
        gsMatrix<> LT; // Lagrange Term
        LT.setZero(n_flat, n_flat);

        gsVector<> zeta = m_dJC.massMatrixSolve(m_lambda); // M^(-1) * m_lambda, some kind of adjoint variable

        index_t iter = 0;
        for(index_t p = 0; p < m_mp->nBoxes(); p++){
            for(index_t i = 0; i < m_dJC.m_detJacBasis.size(p); i++){
                // gsDebugVar(m_dJC.hessD(p,i).rows());
                // gsDebugVar(m_dJC.hessD(p,i).cols());
                LT += m_dJC.hessD(p,i)*zeta[iter++];
            }
        }

        gsMatrix<> tmp = mapMatrix(m_dJC.space_mapper(),LT);
        gsMatrix<> LTTagged = mapMatrixToTagged(m_dJC.space_mapper(),tmp);
        LT = mapMatrix(m_dJC.space_mapper(),tmp);

        // Compute hessian of objective
        gsMatrix<> hessO = hessObj(hessObjTagged);

        // Add Lagrange terms
        hessObjTagged += LTTagged;

        return hessO + LT;
    } else { // If there are no constraints
        return hessObj(hessObjTagged);
    }
}
