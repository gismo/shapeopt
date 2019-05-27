#include <gismo.h>
#include "gsShapeOptProblem.h"

// Implement
// FIXIT: Implement functionality?
gsShapeOptProblem::gsShapeOptProblem(gsMultiPatch<>* mp, gsShapeOptLog* slog):
    m_mp(mp),
    m_dJC(mp),
    m_mappers(mp->targetDim()),
    m_log(slog)
{
    gsInfo << "\n Constructor 1 \n";
};

// Implement
gsShapeOptProblem::gsShapeOptProblem(gsMultiPatch<>* mp, std::vector< gsDofMapper > mappers, gsShapeOptLog* slog):
    m_mp(mp), m_mappers(mappers), m_dJC(m_mp), m_log(slog)
{
    gsInfo << "\n Constructor 2 \n";
    n_free = m_paramMethod->n_free;
    n_flat = m_paramMethod->n_flat;
    n_tagged = m_paramMethod->n_tagged;
    n_cps = m_paramMethod->n_cps;
    setupOptParameters();
};

void gsShapeOptProblem::setupOptParameters()
{
    // The desing variables is the free variables
    m_curDesign = getTagged();
    m_numDesignVars = n_tagged;

    // Default behaviour is without design bounds
    m_desLowerBounds.setConstant(n_tagged, -1e9);
    m_desUpperBounds.setConstant(n_tagged, 1e9);

    // Default constraints are the gsDetJacConstraints
    // Call once to setup solver
    m_dJC.evalCon();

    m_numConstraints = m_dJC.numConstraints();
    m_conLowerBounds = m_dJC.getLowerBounds();
    m_conUpperBounds = m_dJC.getUpperBounds();

    // compute jac structure, default is without sparsity
    computeJacStructure();
};

void gsShapeOptProblem::computeJacStructure()
{
    // Default is without sparsity
    gsIpOptSparseMatrix J = jacobDetJac(); // Jacobiant of detJ constraints

    m_numConJacNonZero = J.nnz();
    m_conJacRows = J.rows();
    m_conJacCols = J.cols();

};

gsVector<> gsShapeOptProblem::evalCon() const
{
    return m_dJC.evalCon();
};

real_t gsShapeOptProblem::evalObj( const gsAsConstVector<real_t> & u ) const
{
    // gsInfo << "evalObj\n";
    updateDesignVariables(u);
    return evalObj();
}

void gsShapeOptProblem::gradObj_into ( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result ) const
{
    // gsInfo << "gradObj_into\n";
    updateDesignVariables(u);
    result = gradObj();
}

void gsShapeOptProblem::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "evalCon_into\n" << std::flush;
    updateDesignVariables(u);
    result = evalCon();
    // gsInfo << "max result " << result.maxCoeff() << "\n";
}

void gsShapeOptProblem::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "jacobCon_into\n";
    updateDesignVariables(u);
    gsIpOptSparseMatrix J = jacobDetJac(); // Jacobiant of detJ constraints
    result = J.values();
}

gsVector<> gsShapeOptProblem::getDesignVariables() const
{
    return m_paramMethod->getTagged();
}

void gsShapeOptProblem::updateDesignVariables(gsVector<> u) const
{
    m_paramMethod->update(u);
}

gsMatrix<> gsShapeOptProblem::jacobDesignUpdate() const
{
    gsMatrix<> out(n_free + n_tagged,n_tagged);
    out.block(0,0,n_free,n_tagged) = m_paramMethod->jacobUpdate(getTagged());
    out.block(n_free,0,n_tagged,n_tagged).setIdentity(n_tagged,n_tagged);
    return out;
}

gsMatrix<> gsShapeOptProblem::mapGradient(gsDofMapper mapper_in, gsMatrix<> mat_in) const
{
    // Set shifts of input mapper
    gsVector<> mapper_in_shifts, tagged_shift;
    mapper_in_shifts.setZero(m_mp->targetDim());
    tagged_shift.setZero(m_mp->targetDim());

    // The first tagged is after all the free cps
    tagged_shift[0] = n_free;

    for(index_t d = 1; d < m_mp->targetDim(); d++){
        mapper_in_shifts[d] = mapper_in_shifts[d-1] + mapper_in.freeSize();
        tagged_shift[d] = tagged_shift[d-1] + m_mappers[d].taggedSize();
    }

    // FIXIT: take this information as input instead..
    bool row = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.rows();
    bool col = mapper_in.freeSize() + mapper_in_shifts[m_mp->targetDim()-1] == mat_in.cols();
    gsMatrix<> mat_out;

    if(row){
        mat_out.setZero(n_free + n_tagged,mat_in.cols());
    } else if (col) {
        mat_out.setZero(mat_in.rows(),n_free + n_tagged);
    } else {
        GISMO_ERROR("Wrong input size in mapMatrix..\n");
    }

    // FIXIT: Calculate all gradient using "setInterfaceCont(0)",
    //      and simply compute this as a permumation (e.g. with a permumation matrix?)

    // Map the free
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
                    mat_out.row(ii + m_paramMethod->m_shift_free[d]) += mat_in.row(ii2);
                } else if (col){
                    mat_out.col(ii + m_paramMethod->m_shift_free[d]) += mat_in.col(ii2);
                }
            }
        }
    }

    // Map the tagged
    for(index_t d = 0; d < m_mp->targetDim(); d++){
        // Iterate through free indices
        for (index_t t = 0; t < m_mappers[d].taggedSize(); t++){
            // Get global index
            index_t ii = m_mappers[d].getTagged()[t];
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
                    mat_out.row(t + tagged_shift[d]) += mat_in.row(ii2);
                } else if (col){
                    mat_out.col(t + tagged_shift[d]) += mat_in.col(ii2);
                }
            }
        }
    }

    return mat_out;

}

// FIXIT sparse structure of dDdc is lost.. Recover possibly?
gsIpOptSparseMatrix gsShapeOptProblem::jacobDetJac() const
{
    gsMatrix<> dDdc = mapGradient(m_dJC.space_mapper(), m_dJC.getJacobian().asDense());
    gsMatrix<> dcdx = jacobDesignUpdate();

    gsMatrix<> dDdx = dDdc*dcdx;
    gsIpOptSparseMatrix out(dDdx,-1); // -1 indicates that it is treated as a dense matrix
    return out;
}

// Implement
void gsShapeOptProblem::runOptimization(index_t maxiter){
    *m_log << "N.o. cps: " << n_cps << "\n";
    *m_log << "N.o. free cps: " << n_free<< "\n";
    *m_log << "N.o. tagged cps: " << n_tagged << "\n";
    *m_log << "N.o. flat cps: " << n_flat << "\n\n";

    // gsInfo << "DoFs for analysis: " << m_stateEq.dbasis.size() << "\n";

    counter2 = 0;
    if (m_log->plotDesign()){
        // Plot m_mp with name design_counter1.pvd
        std::string name = "design";
        m_log->plotInParaview(*m_mp,name,counter2);
    }
    // Run optimization
    for (index_t i = 0; i < maxiter; i++){
        counter1 = 0;

        // Update parametrization
        *m_log << " Objective function before updating: " << evalObj() << "\n";
        *m_log << " Min d before updating: " << -m_dJC.evalCon().maxCoeff() << "\n\n";
        m_paramMethod->updateAndReset();
        m_paramMethod->computeMap();

        // FIXIT: log whether param was succesful
        *m_log << " Objective function after updating: " << evalObj() << "\n";
        *m_log << " Min d after updating: " << -m_dJC.evalCon().maxCoeff() << "\n\n";

        std::string nameAU = "cps_afterUpdate";
        m_log->saveVec(getFlat(),nameAU ,counter2);

        gsMultiPatch<> dJ = m_dJC.getDetJ();
        std::string namedj = "detJ";
        m_log->plotMultiPatchOnGeometry(*m_mp,dJ,namedj,counter2);

        namedj = "detJ_act";
        m_log->plotActiveMultiPatchOnGeometry(*m_mp,dJ,m_dJC.eps(),namedj,counter2);

        // Solve the current optimization problem
        solve();

        if (m_log->plotDesign()){
            // Plot m_mp with name design_counter1.pvd
            std::string name = "design";
            m_log->plotInParaview(*m_mp,name,counter2);
        }

        counter2++;

        // Check if parametrization is good (larger than double of m_eps)
        // real_t minD = -m_dJC.evalCon().maxCoeff();
        // if (minD > 2*m_dJC.m_eps){
            // gsInfo << "\n\nFinal Solution is found, with min d of "
            // << minD << "\n\n";
            // return;
        // }
    }
}

bool gsShapeOptProblem::intermediateCallback() {
    if (m_log->saveCps()){
        std::string name = "cps";
        m_log->saveVec(getFlat(),name,counter2,counter1++);
        m_log->logObj(evalObj()); // Recalculated, is there a better way?
    }
    return true;
}
