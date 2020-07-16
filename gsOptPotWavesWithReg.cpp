#include <gismo.h>
#include "gsOptPotWavesWithReg.h"


gsOptPotWavesWithReg::gsOptPotWavesWithReg(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptProblem> sopt, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, real_t quA, index_t quB, real_t eps, bool useConstraints, bool glueInterfaces, bool usePow):
    gsShapeOptWithReg(mp, sopt, numRefine, slog, quA, quB, eps, glueInterfaces, usePow),
    m_useConstraints(useConstraints)
{
    // setupMappers called in gsShapeOptWithReg constructor 
    
    m_constraint = memory::make_shared( new gs2NormConstraints(m_mp, m_mappers, m_a, m_b) );

    setupOptParameters();
    updateDesignBounds();

}

// Overload to get 2NormConstraints
void gsOptPotWavesWithReg::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const 
{
    //gsInfo << "EVALCON\n" << std::flush;
    updateDesignVariables(u);
    result = m_constraint->evalCon();

    //gsDebugVar(result);
    //gsDebugVar(m_conLowerBounds);
    //gsDebugVar(m_conUpperBounds);
    //exit(0);


}

void gsOptPotWavesWithReg::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const 
{

    //gsInfo << "JACOBCON\n" << std::flush;
    updateDesignVariables(u);

    gsIpOptSparseMatrix J = m_constraint->getJacobian();

    result = J.values();

}

// Method to setup constraints
// Overloaded to use gs2NormConstraints
void gsOptPotWavesWithReg::setupConstraints()
{
    if (m_useConstraints)
    {
        m_numConstraints = m_constraint->numConstraints();
        m_conLowerBounds = m_constraint->getLowerBounds();
        m_conUpperBounds = m_constraint->getUpperBounds();

        //compute jac structure
        computeJacStructure();
    }
    else 
    {
        m_numConstraints = 0;
    }

}

void gsOptPotWavesWithReg::updateDesignBounds()
{
    for (index_t i = 0; i < n_free; i++)
    {
        index_t d=-1;
        if( i < m_winslow->m_shift_free[1])     // Direction x
        {
            continue;
        }
        else if( i < m_winslow->m_shift_free[2])// Direction y
        {
            continue;
        }
        else                            // Direction z
        {
            // In the z direction we just have and upper and lower bound.
            m_desUpperBounds[i] = 0; 

        }
    }

}
    
void gsOptPotWavesWithReg::computeJacStructure()
{
    // Default is without sparsity
    gsIpOptSparseMatrix J = m_constraint->getJacobian(); // Jacobiant of detJ constraints

    m_numConJacNonZero = J.nnz();
    m_conJacRows = J.rows();
    m_conJacCols = J.cols();

};

bool gsOptPotWavesWithReg::intermediateCallback() {
    // FIXIT objective evaluated again!
    real_t obj = m_opt->evalObj();
    real_t winslow = m_winslow->evalObj();
    // real_t gradn = gradObj().norm();

    gsVector<> v(2);
    v << obj, winslow;

    if (m_log->saveCps()){
        std::string name = "cps";
        m_log->saveVec(m_winslow->getFlat(),name,counter2,counter1);
        m_log->logObj(v); // Recalculated, is there a better way?

        gsMultiPatch<> ur = m_opt->getUR();
        gsMultiPatch<> ui = m_opt->getUI();

        if (counter1 % 5 == 0)
        {
            name = "paraview/mp";
            m_log->plotInParaview(*m_mp,name,counter1);

            name = "paraview/ur";
            m_log->plotMultiPatchOnGeometry(*m_mp,ur,name,counter1);

            name = "paraview/ui";
            m_log->plotMultiPatchOnGeometry(*m_mp,ui,name,counter1);

            name = "xml/ur";
            m_log->saveAsXML(ur,name,counter1);

            name = "xml/ui";
            m_log->saveAsXML(ui,name,counter1);
        }

        counter1++;
    }

    return true;
}

void gsOptPotWavesWithReg::runOptimization()
{
    *m_log << "eps = " << m_eps << "\n";
    std::string name = "paraview/initial_design";
    m_log->plotInParaview(*m_mp,name);

    solve();

    name = "paraview/final_design";
    m_log->plotInParaview(*m_mp,name);

}
    
    
