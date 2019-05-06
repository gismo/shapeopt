#include <gismo.h>
#include "gsOptParamMethod.h"
using namespace gismo;

// FIXIT: default argument does not work
gsOptParamMethod::gsOptParamMethod(gsMultiPatch<>* mpin, bool use_dJC = true):
    gsParamMethod(mpin), m_dJC(mpin), use_detJacConstraint(use_dJC)
{
    setupOptParameters();
};

gsOptParamMethod::gsOptParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers):
        gsParamMethod(mpin, mappers), m_dJC(mpin)
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
        m_conLowerBounds = m_dJC.getLowerBounds();
        m_conUpperBounds = m_dJC.getUpperBounds();
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

// Implement
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

// Implement
void gsOptParamMethod::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "jacobCon_into\n" << std::flush;
    updateFree(u);
    if (use_detJacConstraint){
        gsIpOptSparseMatrix J = mapMatrix(m_dJC.space_mapper(),m_dJC.getJacobian());
        result = J.values();
    } else {
        return;
    }

};

// Implement
void gsOptParamMethod::computeJacStructure()
{
    // Use the sparsity provided by gsDetJacConstraint m_dJC
    if (use_detJacConstraint){
        gsIpOptSparseMatrix J = mapMatrix(m_dJC.space_mapper(),m_dJC.getJacobian());
        gsDebugVar(J.nrows());
        gsDebugVar(J.ncols());

        m_numConJacNonZero = J.nnz();
        m_conJacRows = J.rows();
        m_conJacCols = J.cols();
    } else {
        m_numConJacNonZero = 0;
    }

};

// FIXIT: looses sparse structure, implement to keep sparsity structure, the
// structure can be predicted by adding two matrices with 1 at nonzeros
gsIpOptSparseMatrix gsOptParamMethod::mapMatrix(gsDofMapper mapper_in, gsIpOptSparseMatrix M) const
{
    gsMatrix<> mat_in = M.asDense();

    gsMatrix<> mat_out = mapMatrix(mapper_in,mat_in);

    gsIpOptSparseMatrix out(mat_out,-1);

    return out;

};

gsMatrix<> gsOptParamMethod::mapMatrix(gsDofMapper mapper_in, gsMatrix<> mat_in) const
{
    // Set shifts of input mapper
    gsVector<> mapper_in_shifts;
    mapper_in_shifts.setZero(m_mp->targetDim());

    for(index_t d = 1; d < m_mp->targetDim(); d++){
        mapper_in_shifts[d] = mapper_in_shifts[d-1] + mapper_in.freeSize();
    }

    // FIXIT: take this information as input instead..
    bool row = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.rows();
    bool col = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.cols();
    gsMatrix<> mat_out;

    if(row){
        mat_out.setZero(n_free,mat_in.cols());
    } else if (col) {
        mat_out.setZero(mat_in.rows(),n_free);
    }

    for(index_t d = 0; d < m_mp->targetDim(); d++){
        // Iterate through free indices
        for (index_t ii = 0; ii < m_mappers[d].freeSize(); ii++){
            // Get a local index
            std::vector<std::pair<index_t,index_t> > result;
            m_mappers[d].preImage(ii, result);

            for(std::vector<std::pair<index_t,index_t>>::iterator it=result.begin(); it != result.end(); ++it)
            {
                // Get local index and patch
                index_t p = it->first;
                index_t i = it->second;

                // Convert to global to find the right column
                index_t ii2 = mapper_in.index(i,p) + mapper_in_shifts[d];

                if (row){ // If right no colms
                    mat_out.row(ii + m_shift_free[d]) += mat_in.row(ii2);
                } else if (col){
                    mat_out.col(ii + m_shift_free[d]) += mat_in.col(ii2);
                }
            }
        }
    }
    return mat_out;

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

void gsOptParamMethod::print(){
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
