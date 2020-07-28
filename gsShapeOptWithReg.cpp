#include <gismo.h>
#include <math.h>
#include "gsShapeOptWithReg.h"
#include "gsOptParam.h"
#include "gsWinslowPow.h"
#include "gsWinslowG0Pow.h"

// Implement
gsShapeOptWithReg::gsShapeOptWithReg(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptProblem> sopt, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, real_t quA, index_t quB,real_t eps, bool glueInterfaces, bool usePow, bool useG0):
    m_glueInterfaces(glueInterfaces),
    m_eps(eps),
    m_mp(mp),
    m_log(slog),
    m_mappers(mp->targetDim()),
    m_opt(sopt)
    // m_optA(mp,numRefine,slog,0,quA,quB,false)
{
    setupMappers();

    if (usePow)
    {
        *m_log << "use Pow!\n";
        gsInfo << "use Pow!\n";
        if (useG0)
        {
            *m_log << "use G0 !\n";
            gsInfo << "use G0 !\n";
            m_winslow = memory::make_shared( new gsWinslowG0Pow(m_mp,m_mappers,false,false,true,0) );
            m_winslow->setMp0(*mp);
        }
        else
            m_winslow = memory::make_shared( new gsWinslowPow(m_mp,m_mappers,false,false,true,0) );
    }
    else
    {
        *m_log << "We do not use Pow!\n";
        m_winslow = memory::make_shared(new gsWinslow(m_mp,m_mappers,false,false,true,0));
    }
    m_winslow->setQuad(quA,quB);


    setupOptParameters();
}

real_t gsShapeOptWithReg::evalObj() const
{
    //gsInfo << "EVALOBJ\n" << std::flush;
    real_t winslow = m_winslow->evalObj();
    //gsInfo << "winslow : " << winslow << "\n" <<  std::flush;
    if ( std::isinf(winslow) )
        return winslow;

    //gsInfo << "SE : quA, quB: " << m_opt->getQuA() << ", " << m_opt->getQuB() << "\n";
    //gsInfo << "Win: quA, quB: " << m_winslow->m_quA << ", " << m_winslow->m_quB << "\n";

    return m_opt->evalObj() + m_eps*m_winslow->evalObj();
}

gsVector<> gsShapeOptWithReg::gradObj() const
{
    //gsInfo << "GRADOBJ\n" << std::flush;

    gsVector<> gradAll = m_opt->gradAll() + m_eps * m_winslow->gradAll();
    gsVector<> out = m_winslow->mapMatrix(m_opt->mapper_grad(),gradAll);
    return out;
};

void gsShapeOptWithReg::setupMappers()
{
    // We remake the mapper from gsShapeOptProblem m_opt but we don't fixed the tagged cps

    gsMultiBasis<> geoBasis(*m_mp);

    if (m_glueInterfaces)
    {
        for (index_t d = 0; d < m_mp->targetDim(); d++)
            geoBasis.getMapper(iFace::glue,m_mappers[d],false); // False means that we do not finalize
    } 
    else 
    {
        for (index_t d = 0; d < m_mp->targetDim(); d++)
        {
            gsDofMapper tmp(geoBasis);
            m_mappers[d] = tmp;
            m_mappers[d].finalize();
        }
        return;
    }

    m_shift_free.setZero(m_mp->targetDim());
    for (index_t d = 0; d < m_mp->targetDim(); d++) // For each dimension
    {
        for (index_t k = 0; k < m_mp->nBoxes(); k++) // For each patch
        {
            for (index_t i = 0; i < m_mp->patch(k).coefsSize(); i++) // For each cps
            {
                // If i,k is neither free nor tagged
                if (!m_opt->mappers()[d].is_free(i,k) and !m_opt->mappers()[d].is_tagged(i,k))
                    m_mappers[d].eliminateDof(i,k);

            }
        }
        m_mappers[d].finalize();

        // Tag again
        for (index_t t = 0; t < m_opt->mappers()[d].taggedSize(); t++)
        {
            index_t ii = m_opt->mappers()[d].getTagged()[t];
            std::vector< std::pair<index_t,index_t> > result;
            m_opt->mappers()[d].preImage(ii, result);

            m_mappers[d].markTagged(result[0].second, result[0].first);
        }

        // DEBUG: The order is not kept. Instead iterate over tagg
        // Tag tagged again
        //for (index_t k = 0; k < m_mp->nBoxes(); k++) // For each patch
        //{
        //    for (index_t i = 0; i < m_mp->patch(k).coefsSize(); i++) // For each cps
        //    {
        //        if (m_opt->mappers()[d].is_tagged(i,k))
        //            m_mappers[d].markTagged(i,k);

        //    }
        //}

        gsInfo << "\nDimension " << d << "\n";
        gsInfo << "size " << m_opt->mappers()[d].size() << ", " << m_mappers[d].size() << "\n";
        gsInfo << "taggedSize " << m_opt->mappers()[d].taggedSize() << ", " << m_mappers[d].taggedSize() << "\n";
        gsInfo << "freeSize " << m_opt->mappers()[d].freeSize() << ", " << m_mappers[d].freeSize() << "\n";
        if (d > 0)
            m_shift_free[d] = m_shift_free[d-1] + m_mappers[d-1].freeSize();


    }

    // OLD CODE FOR ANTENNA
    // // Get mappers from multibasis with interfaces glued
    // gsMultiBasis<> geoBasis(*m_mp);
    // geoBasis.getMapper(iFace::glue,m_mappers[0],false); // False means that we do not finalize
    // geoBasis.getMapper(iFace::glue,m_mappers[1],false); // False means that we do not finalize
    //
    // // Fix y coordinate for all boundaries
    // for (index_t i = 0; i < m_mp->nBoundary(); i ++){
    //     patchSide ps = m_mp->boundaries()[i];
    //     gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
    //     // gsInfo << ps << " is fixed in y direction \n\n";
    //     // gsInfo << "boundaryDofs on patch " << ps.patch << "\n"<< boundaryDofs << "\n\n";
    //     m_mappers[1].markBoundary(ps.patch,boundaryDofs);
    //
    //     if (m_mp->nBoxes() > 3){
    //         if (ps.patch == 4){
    //             m_mappers[0].markBoundary(ps.patch,boundaryDofs);
    //         }
    //     } else {
    //         // If there is less than 3 patches then all boundaries is fixed
    //         m_mappers[0].markBoundary(ps.patch,boundaryDofs);
    //     }
    // }
    //
    // // Fix boundary of the fixedPatch
    // gsMultiPatch<> fixedGeom(m_mp->patch(fixedPatch));
    //
    // // Finalize mappers
    // m_mappers[0].finalize();
    // m_mappers[1].finalize();
    //
    // gsInfo << "n controlpoints: " << m_mappers[0].mapSize() << "\n";
    //
    // // Tag fixed bnds
    // for(index_t i = 0; i < fixedGeom.nBoundary(); i++){
    //     patchSide ps = fixedGeom.boundaries()[i];
    //     gsVector<unsigned> boundaryDofs = m_mp->basis(fixedPatch).boundary(ps);
    //     for (index_t j = 0; j < boundaryDofs.size(); j ++){
    //         m_mappers[0].markTagged(boundaryDofs[j],fixedPatch); // Mark xcoord of boundaries of fixed patch
    //         m_mappers[1].markTagged(boundaryDofs[j],fixedPatch); // Mark ycoord of boundaries of
    //
    //     }
    // }
    //
    // // count number of free dofs
    // n_free = 0;
    // for(index_t i = 0; i < m_mp->targetDim(); i++){
    //     n_free += m_mappers[i].freeSize();
    // }
}

real_t gsShapeOptWithReg::evalObj( const gsAsConstVector<real_t> & u ) const
{
    updateDesignVariables(u);
    return evalObj();
}

void gsShapeOptWithReg::gradObj_into ( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result ) const
{
    result = gradObj();
}

void gsShapeOptWithReg::updateDesignVariables(gsVector<> u) const
{
    m_winslow->updateFree(u);
}

void gsShapeOptWithReg::setupOptParameters()
{
    n_free = m_winslow->n_free;
    n_flat = m_winslow->n_flat;
    n_tagged = m_winslow->n_tagged;
    n_cps = m_winslow->n_cps;

    // The desing variables is the free variables
    m_curDesign = m_winslow->getFree();
    m_numDesignVars = n_free;

    // Default behaviour is without design bounds
    m_desLowerBounds.setConstant(n_free, -1e9);
    m_desUpperBounds.setConstant(n_free, 1e9);

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

void gsShapeOptWithReg::setupConstraints()
{
    m_numConstraints = 0;
    m_numConJacNonZero = 0;
}

void gsShapeOptWithReg::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    GISMO_ERROR("evalCon called");
}

void gsShapeOptWithReg::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    GISMO_ERROR("jacobCon called");
}

bool gsShapeOptWithReg::intermediateCallback() {
    real_t obj = m_opt->evalObj();
    real_t winslow = m_winslow->evalObj();
    // real_t gradn = gradObj().norm();

    gsVector<> v(2);
    // gsInfo << "norm of grad: " << gradObj().norm() << "\n";
    v << obj, winslow;

    // if (counter1 % 10 == 0)
    // {
    //     std::string name = "geo";
    //     m_log->plotInParaview(*m_mp,name,counter1);
    // }

    if (m_log->saveCps()){
        std::string name = "cps";
        m_log->saveVec(m_winslow->getFlat(),name,counter2,counter1++);
        m_log->logObj(v); // Recalculated, is there a better way?

    }
    return true;
}

void gsShapeOptWithReg::runOptimization(index_t maxiter)
{
    *m_log << "N.o. cps: " << n_cps << "\n";
    *m_log << "N.o. free cps: " << n_free<< "\n";
    *m_log << "N.o. tagged cps: " << n_tagged << "\n";
    *m_log << "N.o. flat cps: " << n_flat << "\n\n";

    *m_log << "m_eps = " << m_eps << "\n\n";

    // gsInfo << "DoFs for analysis: " << m_stateEq.dbasis.size() << "\n";

    counter2 = 0;

    real_t k = 1.0/sqrt(2); // factor to decrease eps with every iteration
    *m_log << "factor to decrease m_eps; k = " << k << "\n\n";

    // Run optimization
    for (index_t i = 0; i < maxiter; i++){
        counter1 = 0;

        *m_log << "m_eps = " << m_eps << "\n";

        // Solve the current optimization problem
        solve();

        counter2++;

        m_eps *= k;
    }
};

void gsShapeOptWithReg::setWinslowQuad(real_t quA, index_t quB)
{
    m_winslow->setQuad(quA, quB);
}

void gsShapeOptWithReg::runOptimizationUntilPosDetJ(index_t maxiter, real_t k, index_t maxRef)
{
    *m_log << "N.o. cps: " << n_cps << "\n";
    *m_log << "N.o. free cps: " << n_free<< "\n";
    *m_log << "N.o. tagged cps: " << n_tagged << "\n";
    *m_log << "N.o. flat cps: " << n_flat << "\n\n";

    *m_log << "Start with m_eps = " << m_eps << "\n";

    real_t init_minDJ15 = m_winslow->minDetJInGaussPts(15);  // 15 additional pts (quA,quB+15)
    real_t init_maxDJ15 = m_winslow->maxDetJInGaussPts(15);  // 15 additional pts (quA,quB+15)
    gsInfo << "Smallest value of detJ mp in 15 extra pts \t: " << init_minDJ15 << "\n";
    gsInfo << "Largest value of detJ mp in 15 extra pts \t: " << init_maxDJ15 << "\n";

    // gsInfo << "DoFs for analysis: " << m_stateEq.dbasis.size() << "\n";

    counter2 = 0;

    *m_log << "   and factor to decrease m_eps; k = " << k << "\n\n";

    *m_log << "---- Run Optimization Until Positive Det J ----\n";

    // Run optimization
    for (index_t i = 0; i < maxiter; i++){
        counter1 = 0;

        *m_log << "m_eps = " << m_eps << "\n";
        gsInfo << "m_eps = " << m_eps << "\n";

        // Solve the current optimization problem
        solve();

        // Save the result in a paraview file
        std::string name = "mp";
        m_log->plotInParaview(*m_mp,name,i);

        // Compute the smallest value of detJ in quadrature pts
        real_t minDJ0 = m_winslow->minDetJInGaussPts(0);    // The pts we use in the optimization
        real_t minDJ15 = m_winslow->minDetJInGaussPts(15);  // 15 additional pts (quA,quB+15)
        *m_log << "Smallest value of detJ mp in quad pts \t\t: " << minDJ0 << "\n";
        *m_log << "Smallest value of detJ mp in 15 extra pts \t: " << minDJ15 << "\n";

        // Check wether to increse quad pts
        if ((minDJ0 - minDJ15)/minDJ15 > 1e-3)
        {
            real_t quA = m_winslow->m_quA;
            real_t quB = m_winslow->m_quB;

            real_t quAnew = quA + 1;
            real_t quBnew = quB + 1;
            m_winslow->setQuad(quAnew, quBnew);

            *m_log << "Relative difference in detJ value: " << (minDJ0 - minDJ15)/minDJ15 << "\n";
            *m_log << "Increase (quA,quB) from (" << quA << ", " << quB << ") to (";
            *m_log << quAnew << ", " << quBnew << ")\n";

        }

        // Increase counter and decrease m_eps
        counter2++;
        m_eps *= k;

        // Snap bnd control points
        gsMultiPatch<> snapped  = (std::dynamic_pointer_cast< gsOptParam >(m_opt))->getSnapped();
        gsMultiPatch<>::Ptr snapped_ptr = memory::make_shared_not_owned(&snapped);

        // Save the snapped domain
        name = "snapped";
        m_log->plotInParaview(snapped,name,i);

        // Construct gsDetJacConstraint
        gsInfo << "construct dJC\n";
        gsDetJacConstraint dJC(snapped_ptr, true); // True means that we use tensor product structure
        index_t neededSteps;

        // Construct temporary Winslow method
        gsWinslow tmp(snapped_ptr, false);
    
        gsInfo << "compute minDetJInGauss\n";
        real_t minDetJ0_snapped = tmp.minDetJInGaussPts(0);
        real_t minDetJ15_snapped = tmp.minDetJInGaussPts(15);

        gsInfo << "detJ value of of snapped :" << minDetJ0_snapped << " \n";
        gsInfo << "detJ value of of snapped :" << minDetJ15_snapped << " \n";
        *m_log << "detJ value of of snapped :" << minDetJ0_snapped << " \n";
        *m_log << "detJ value of of snapped :" << minDetJ15_snapped << " \n";

        if (minDetJ15_snapped < 0)
        {
            gsInfo << "Negative detJ of snapped found, continues!\n\n";
            *m_log << "Negative detJ of snapped found, continues!\n\n";
            continue;
        }

        gsInfo << "Positive detJ of snapped found, continues!\n\n";
        gsInfo << "compute detJ coefs\n";
        real_t minD = dJC.provePositivityOfDetJ_TP(neededSteps, maxRef);

        if (minD > 0)
        {
            gsInfo << "Optimization ended with minD = " << minD << "\n";
            *m_log << "Optimization ended with minD = " << minD << "\n";
            return;
        }
    }
};
