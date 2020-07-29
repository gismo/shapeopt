#include <gismo.h>
#include "gsWinslowPow.h"
#include "gsWinslowG0Pow.h"
#include "gsDetJacConstraint.h"
#include "gsOptPotWavesWithReg2nd.h"
#include "gsDetJacConstraint.h"



gsOptPotWavesWithReg2nd::gsOptPotWavesWithReg2nd(gsMultiPatch<>::Ptr mp, gsMultiPatch<>::Ptr center, gsShapeOptProblem::Ptr sopt, index_t numRefine, gsShapeOptLog::Ptr slog, real_t quA, index_t quB, real_t eps,  bool glueInterfaces, bool usePow, bool useG0):
    gsShapeOptWithReg(mp, sopt, numRefine, slog, quA, quB, eps, glueInterfaces, usePow, useG0),
    m_center(center)
{
    // setupMappers called in gsShapeOptWithReg2nd constructor 
    
    // Setup center param method
    if (usePow)
        if (useG0)
        {
            m_center_winslow = memory::make_shared( new gsWinslowG0Pow(m_center,false,false,true,0) );
            m_center_winslow->setMp0(*center);
        }
        else
            m_center_winslow = memory::make_shared( new gsWinslowPow(m_center,false,false,true,0) );
    else 
        m_center_winslow = memory::make_shared( new gsWinslow(m_center,false,false,true,0) );

    m_winslow->setQuad(quA,quB);
    // Compute map from control points to bnd of center
    n_constraints = 0;

    // For each cps in center
    for(index_t j = 0; j < m_center->patch(0).coefsSize(); j++)
    {
        gsVector<> coefc = m_center->patch(0).coef(j);

        bool flag = false;
        // For each cps
        for (index_t p = 0; p < m_mp->nBoxes(); p++)
        {
            for(index_t i = 0; i < m_mp->patch(p).coefsSize(); i++)
            {
                gsVector<> coef = m_mp->patch(p).coef(i);

                if ((coef - coefc).norm() == 0)
                {
                    for ( index_t d = 0; d < 3; d++)
                    {
                        index_t ii = m_mappers[d].index(i,p);

                        if (m_mappers[d].is_free_index(ii))
                        {
                            n_constraints++;
                        }
                    }
                    flag = true;
                    break;
                }
            }
            if (flag)
                break;
        }
    }

    setupOptParameters();
    updateDesignBounds();

    m_constraint_matrix.setZero(n_constraints, m_numDesignVars);
    index_t iflat = 0;
    index_t iconst = 0;

    gsVector<> center_flat = m_center_winslow->getFlat();
    gsVector<> free = m_winslow->getFree();

    for(index_t j = 0; j < m_center->patch(0).coefsSize(); j++)
    {
        gsVector<> coefc = m_center->patch(0).coef(j);

        bool flag = false;
        // For each cps
        for (index_t p = 0; p < m_mp->nBoxes(); p++)
        {
            for(index_t i = 0; i < m_mp->patch(p).coefsSize(); i++)
            {
                gsVector<> coef = m_mp->patch(p).coef(i);

                if ((coef - coefc).norm() == 0)
                {
                    for ( index_t d = 0; d < 3; d++)
                    {
                        index_t ii = m_mappers[d].index(i,p);

                        if (m_mappers[d].is_free_index(ii))
                        {
                            index_t kk = ii + m_winslow->m_shift_free[d];

                            m_constraint_matrix(iconst,kk) = 1;
                            m_constraint_matrix(iconst,n_free + j + m_center_winslow->m_shift_flat[d]) = -1;
                            iconst++;
                        }
                        else
                        {
                            m_desLowerBounds[n_free + j + m_center_winslow->m_shift_flat[d] ] = coef[d];
                            m_desUpperBounds[n_free + j + m_center_winslow->m_shift_flat[d]] = coef[d];
                        }

                    }
                    flag = true;
                    break;
                }
            }
            if (flag)
                break;
        }
    }

    gsIpOptSparseMatrix J(m_constraint_matrix,0);
    m_J = J;

    computeJacStructure();
    
    //gsMatrix<> disp(m_numDesignVars,3);
    //disp << m_desLowerBounds, m_curDesign, m_desUpperBounds;
    //gsDebugVar(disp);
    


}

void gsOptPotWavesWithReg2nd::setupOptParameters()
{
    n_free = m_winslow->n_free;
    gsDebugVar(n_free);
    n_flat = m_winslow->n_flat;
    n_tagged = m_winslow->n_tagged;
    n_cps = m_winslow->n_cps;

    n_flat_center = m_center_winslow->n_flat;

    // The desing variables is the free variables
    m_numDesignVars = n_free + n_flat_center;

    m_curDesign.setZero(m_numDesignVars,1);
    m_curDesign.block(0,0,n_free,1) = m_winslow->getFree();
    m_curDesign.block(n_free,0,n_flat_center,1) = m_center_winslow->getFlat();

    // Default behaviour is without design bounds
    m_desLowerBounds.setConstant(m_numDesignVars, -1e9);
    m_desUpperBounds.setConstant(m_numDesignVars, 1e9);

    // Set shifts of input mapper
    gsVector<> shifts = m_winslow->shift_free();

    index_t iter = 0;

    gsVector<> tagged = m_opt->m_paramMethod->getTagged();
    gsVector<> free = m_winslow->getFree();
    // Chose correct bounds for tagged cps
    for(index_t d = 0; d < m_mp->targetDim(); d++){
        // Iterate through tagged indices
        for (index_t t = 0; t < m_opt->mappers()[d].taggedSize(); t++){
            // Get global index
            index_t ii_mopt = m_opt->mappers()[d].getTagged()[t];

            std::vector< std::pair< index_t, index_t > > result;
            m_opt->mappers()[d].preImage(ii_mopt,result);

            index_t ii = m_mappers[d].index(result[0].second, result[0].first);

            index_t jj = ii + shifts[d];

            if (tagged[iter] - m_opt->lowerDesignVar(iter) < 0)
            {
                gsDebugVar(tagged[iter]);
                gsDebugVar(m_opt->lowerDesignVar(iter));
                GISMO_ERROR("ERROR WITH LOWER DES VAR");
            }
            if (-tagged[iter] + m_opt->upperDesignVar(iter) < 0)
            {
                gsDebugVar(tagged[iter]);
                gsDebugVar(m_opt->upperDesignVar(iter));
                GISMO_ERROR("ERROR WITH UPPER DES VAR");
            }

            m_desLowerBounds[jj] = m_opt->lowerDesignVar(iter);
            m_desUpperBounds[jj] = m_opt->upperDesignVar(iter);

            iter++;
        }
    }
     //gsMatrix<> out(m_numDesignVars,2);
     //gsVector<> free = m_winslow->getFree();
     //out << free - m_desLowerBounds, m_desUpperBounds - free;
     //gsInfo << out;

    setupConstraints();
}

void gsOptPotWavesWithReg2nd::setupConstraints()
{
    m_numConstraints = n_constraints;
    m_conLowerBounds.setZero(n_constraints);
    m_conUpperBounds.setZero(n_constraints);

    //compute jac structure
    computeJacStructure();

}

void gsOptPotWavesWithReg2nd::updateDesignBounds()
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

    for (index_t i = 0; i < n_flat_center; i++)
    {
        index_t d=-1;
        if( i < m_center_winslow->m_shift_flat[1])     // Direction x
        {
            //gsInfo << "Skip dir 0 : " << i << "\n";
            continue;
        }
        else if( i < m_center_winslow->m_shift_flat[2])// Direction y
        {
            //gsInfo << "Skip dir 1 : " << i << "\n";
            continue;
        }
        else                            // Direction z
        {
            // In the z direction we just have and upper and lower bound.
            m_desUpperBounds[n_free + i] = 0; 

        }
    }

}

bool gsOptPotWavesWithReg2nd::intermediateCallback() 
{
    // FIXIT objective evaluated again!
    real_t obj = m_opt->evalObj();
    //real_t obj = 0;
    real_t winslow = m_winslow->evalObj() + m_center_winslow->evalObj();
    // real_t gradn = gradObj().norm();

    gsVector<> v(2);
    v << obj, winslow;

    if (m_log->saveCps()){
        std::string name = "cps";
        m_log->saveVec(m_winslow->getFlat(),name,counter2,counter1);
        m_log->logObj(v); // Recalculated, is there a better way?

        gsMultiPatch<> ur = m_opt->getUR();
        gsMultiPatch<> ui = m_opt->getUI();

        if (counter1 % 50 == 0)
        {
            name = "paraview/mp";
            m_log->plotInParaview(*m_mp,name,counter1);

            name = "paraview/center";
            m_log->plotInParaview(*m_center,name,counter1);

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

void gsOptPotWavesWithReg2nd::runOptimization()
{
    *m_log << "eps = " << m_eps << "\n";
    updateDesignVariables(m_curDesign);

    std::string name = "paraview/initial_design";
    m_log->plotInParaview(*m_mp,name);

    name = "paraview/initial_center";
    m_log->plotInParaview(*m_center,name);

    solve();

    *m_log << "min detJ in gauss pts: " << m_winslow->minDetJInGaussPts() << "\n";
    *m_log << "min detJ in gauss pts center: " << m_center_winslow->minDetJInGaussPts() << "\n";

    *m_log << "min detJ in more gauss pts: " << m_winslow->minDetJInGaussPts(15) << "\n";
    *m_log << "min detJ in more gauss pts center: " << m_center_winslow->minDetJInGaussPts(15) << "\n";

    name = "paraview/final_design";
    m_log->plotInParaview(*m_mp,name);

    gsDetJacConstraint dJC(m_mp,true);
    name = "paraview/detJ";
    m_log->plotDetJ(dJC,name);
    

}

void gsOptPotWavesWithReg2nd::computeJacStructure()
{
    m_numConJacNonZero = m_J.nnz();
    m_conJacRows = m_J.rows();
    m_conJacCols = m_J.cols();

};
    
void gsOptPotWavesWithReg2nd::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const 
{
    //gsInfo << "EVALCON\n" << std::flush;
    result = m_constraint_matrix * u;

    //gsDebugVar(result);
    //gsDebugVar(m_conLowerBounds);
    //gsDebugVar(m_conUpperBounds);
    //exit(0);


}

void gsOptPotWavesWithReg2nd::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const 
{
    result = m_J.values();
}

real_t gsOptPotWavesWithReg2nd::evalObj() const
{
    //gsInfo << "EVALOBJ\n" << std::flush;
    real_t winslow = m_winslow->evalObj() ;
    if ( std::isinf(winslow) )
        return winslow;

    winslow += m_center_winslow->evalObj();
    if ( std::isinf(winslow) )
        return winslow;

    //gsInfo << "SE : quA, quB: " << m_opt->getQuA() << ", " << m_opt->getQuB() << "\n";
    //gsInfo << "Win: quA, quB: " << m_winslow->m_quA << ", " << m_winslow->m_quB << "\n";

    return m_opt->evalObj() + m_eps*( winslow );
    //return m_eps*( winslow );
}

gsVector<> gsOptPotWavesWithReg2nd::gradObj() const
{
    //gsInfo << "GRADOBJ\n" << std::flush;
    gsVector<> out(m_numDesignVars);

    gsVector<> gradAll = m_opt->gradAll() + m_eps * m_winslow->gradAll();
    //gsVector<> gradAll = m_eps * m_winslow->gradAll();

    out.segment(0,n_free) = m_winslow->mapMatrix(m_opt->mapper_grad(),gradAll);
    out.segment(n_free,n_flat_center) = m_eps*m_center_winslow->gradAll();

    return out;

};

void gsOptPotWavesWithReg2nd::updateDesignVariables(gsVector<> u) const
{
    m_winslow->updateFree(u.segment(0,n_free));
    m_center_winslow->updateFlat(u.segment(n_free,n_flat_center));
}

void gsOptPotWavesWithReg2nd::setWinslowQuad(real_t quA, index_t quB)
{
    m_winslow->setQuad(quA, quB);
    m_center_winslow->setQuad(quA, quB);
}
