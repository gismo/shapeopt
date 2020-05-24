#include <gismo.h>
#include "gsShapeOptProblem.h"
#include "gsWinslow.h"
#include "gsAffineOptParamMethod.h"

// FIXIT: Implement functionality?
gsShapeOptProblem::gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptLog> slog, bool useDetJCons):
    m_mp(mp),
    m_mappers(mp->targetDim()),
    m_log(slog),
    m_useDetJCons(useDetJCons)
{
    gsInfo << "\n Constructor 1 \n";
    if (m_useDetJCons)
        m_dJC = memory::make_shared(new gsDetJacConstraint(mp));

};

gsShapeOptProblem::gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, memory::shared_ptr<gsShapeOptLog> slog, memory::shared_ptr<gsConstraint> constraint):
    m_mp(mp),
    m_mappers(mp->targetDim()),
    m_log(slog),
    m_dJC(constraint)
{
    gsInfo << "\n Constructor 1 \n";
};

gsShapeOptProblem::gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, std::vector< gsDofMapper > mappers, memory::shared_ptr<gsShapeOptLog> slog):
    m_mp(mp), m_mappers(mappers), m_log(slog)
{
    gsInfo << "\n Constructor 2 \n";
    m_dJC = memory::make_shared(new gsDetJacConstraint(mp));
};

gsShapeOptProblem::gsShapeOptProblem(memory::shared_ptr<gsMultiPatch<>> mp, std::vector< gsDofMapper > mappers, memory::shared_ptr<gsShapeOptLog> slog, memory::shared_ptr<gsConstraint> constraint):
    m_mp(mp), m_mappers(mappers), m_dJC(constraint), m_log(slog)
{
    m_useOwnCons = true;
    gsInfo << "\n Constructor 3 \n";
    setupOptParameters();
};

void gsShapeOptProblem::setupOptParameters()
{
    n_free = m_paramMethod->n_free;
    n_flat = m_paramMethod->n_flat;
    n_tagged = m_paramMethod->n_tagged;
    n_cps = m_paramMethod->n_cps;

    // The desing variables is the free variables
    m_curDesign = getTagged();
    m_numDesignVars = n_tagged;

    // Default behaviour is without design bounds
    setupDesignBounds();
    // m_desLowerBounds.setConstant(n_tagged, -1e9);
    // m_desUpperBounds.setConstant(n_tagged, 1e9);

    setupConstraints();
};

void gsShapeOptProblem::setupConstraints()
{
    // Default constraints are the gsDetJacConstraints
    // Call once to setup solver
    if (m_useDetJCons) {
        m_dJC->evalCon();

        m_numConstraints = m_dJC->numConstraints();
        m_conLowerBounds = m_dJC->getLowerBounds();
        m_conUpperBounds = m_dJC->getUpperBounds();

        //compute jac structure, default is without sparsity
        computeJacStructure();
    } else {
        m_numConstraints = 0;
        m_numConJacNonZero = 0;

    }

}

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
    return m_dJC->evalCon();
};

real_t gsShapeOptProblem::evalObj( const gsAsConstVector<real_t> & u ) const
{
    // gsInfo << "evalObj\n";
    if (!updateDesignVariables(u))
    {
        gsInfo << "return inf\n";
        return std::numeric_limits<double>::infinity();
    }
    return evalObj();
}

void gsShapeOptProblem::gradObj_into ( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result ) const
{
    // gsInfo << " desVars: " << m_numDesignVars << "\n";
    // gsInfo << "gradObj, u = " << u << "\n" << std::flush;
    if (!updateDesignVariables(u)){
        // GISMO_ERROR("UPDATE WENT WRONG IN GRADOBJINTO\n");
        gsInfo << "\n ============================= \n";
        gsInfo << "GRAD OBJ FAILED\n";
        gsInfo << "============================= \n";
    }
    // gsInfo << "gradObj_into\n" << std::flush;
    // gsInfo << "size of gradObj" << gradObj().rows() << "\n\n" << std::flush;
    result = gradObj();
}

void gsShapeOptProblem::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "evalCon_into\n" << std::flush;
    // gsInfo << "evalcon, u = " << u << "\n" << std::flush;
    if (!updateDesignVariables(u)){
        // GISMO_ERROR("UPDATE WENT WRONG IN EVALCONINTO\n");
        gsVector<> tagged = m_paramMethod->getTagged();
        m_paramMethod->updateTagged(u);
        result = evalCon();
        m_paramMethod->updateTagged(tagged);
        gsInfo << "\n ============================= \n";
        gsInfo << "EVAL CON FAILED\n";
        gsInfo << "============================= \n";
    }
    result = evalCon();
}

void gsShapeOptProblem::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "jacobCon_into\n";
    if (!updateDesignVariables(u)){
        // GISMO_ERROR("UPDATE WENT WRONG IN JACOBCONINTO\n");
        gsInfo << "\n ============================= \n";
        gsInfo << "JACOB CON FAILED\n";
        gsInfo << "============================= \n";
    }
    gsIpOptSparseMatrix J = jacobDetJac(); // Jacobiant of detJ constraints
    result = J.values();
}

gsVector<> gsShapeOptProblem::getDesignVariables() const
{
    return m_paramMethod->getTagged();
}

// Returns false if the update failed
bool gsShapeOptProblem::updateDesignVariables(gsVector<> u) const
{
    return m_paramMethod->update(u);
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
    gsMatrix<> dDdc = mapGradient(m_dJC->space_mapper(), m_dJC->getJacobian().asDense());
    gsMatrix<> dcdx = jacobDesignUpdate();

    gsMatrix<> dDdx = dDdc*dcdx;
    gsIpOptSparseMatrix out(dDdx,-1); // -1 indicates that it is treated as a dense matrix
    return out;
}

void gsShapeOptProblem::runOptimization(index_t maxiter, bool uniRef)
{

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
        *m_log << " Min d before updating: " << m_dJC->evalCon().minCoeff() << "\n\n";
        gsInfo << " Min d before updating: " << m_dJC->evalCon().minCoeff() << "\n\n";

        gsVector<> flat = m_paramMethod->getFlat();

        //gsMultiBasis<> intBas(*m_mp);
        // intBas.uniformRefine();
        //(std::dynamic_pointer_cast< gsAffineOptParamMethod >(m_paramMethod))->setIntegrationBasis(intBas);


        bool status = m_paramMethod->updateAndReset();
        // bool status = true;

        // m_paramMethod->refineBasedOnDetJ(0,dynamic_cast< gsDetJacConstraint*>(m_dJC));
        m_dJC->setup();
        real_t mind = m_dJC->evalCon().minCoeff();
        *m_log << " Min d after updating: " << mind << "\n\n";
        gsInfo << " Min d after updating: " << mind << "\n\n";

        gsMultiPatch<> dJ = (std::dynamic_pointer_cast< gsDetJacConstraint >(m_dJC))->getDetJSurface();
        std::string namedj = "detJSurf";
        m_log->plotInParaview(dJ,namedj,counter2);

		if (uniRef)
        	m_dJC->refineUntilPositiveUniformly(7,0);
		else 
        	m_dJC->refineUntilPositive(7,0);
				
        // m_paramMethod->updateAndReset();
        // gsMultiPatch<> dJ = (dynamic_cast< gsDetJacConstraint* >(m_dJC))->getDetJ();
        // std::string namedj = "detJ";
        // m_log->plotMultiPatchOnGeometry(*m_mp,dJ,namedj,ct+10*counter2);
        mind = m_dJC->evalCon().minCoeff();
        *m_log << " Min d after refinement: " << mind << "\n";
        gsInfo << " Min d after refinement: " << mind << "\n";
        // m_dJC->refineUntilPositive(1,mind*1.01);
        // mind = m_dJC->evalCon().minCoeff();
        // *m_log << " Min d after 2nd refinement: " << mind << "\n";
        // gsInfo << " Min d after 2nd refinement: " << mind << "\n";

        m_dJC->setEps(0.25*mind);

        gsMultiBasis<> bas = (std::dynamic_pointer_cast< gsDetJacConstraint >(m_dJC))->m_detJacBasis;
        namedj = "detJBasis";
        m_log->plotMultiBasisOnGeometry(*m_mp,bas,namedj,counter2);

        // m_dJC->refineUntilPositive(5);
        setupOptParameters();
        print();

        // *m_log << "Refined m_dJC <5 times, number of constraints are now " << m_numConstraints << "\n";

        if (not status){
            // GISMO_ERROR("m_paramMethod failed");
        }

        m_paramMethod->computeMap();

        // FIXIT: log whether param was succesful
        *m_log << " Objective function after updating: " << evalObj() << "\n";
        *m_log << " Min d after updating: " << m_dJC->evalCon().minCoeff() << "\n\n";

        m_paramMethod->update(m_curDesign);
        *m_log << " Min d after updating with affine: " << m_dJC->evalCon().minCoeff() << "\n\n";
        gsInfo << " Min d after updating with affine: " << m_dJC->evalCon().minCoeff() << "\n\n";

        std::string nameAU = "cps_afterUpdate";
        m_log->saveVec(getFlat(),nameAU ,counter2);

        // gsMultiPatch<> dJ = m_dJC.getDetJ();
        // std::string namedj = "detJ";
        // m_log->plotMultiPatchOnGeometry(*m_mp,dJ,namedj,counter2);

        // namedj = "detJ_act";
        // m_log->plotActiveMultiPatchOnGeometry(*m_mp,dJ,m_dJC.eps(),namedj,counter2);

        // Solve the current optimization problem
        solve();

        *m_log << " Max Lagrange multipliers from shapeopt " << m_lambda.maxCoeff() << "\n";
        *m_log << " Min Lagrange multipliers from shapeopt " << m_lambda.minCoeff() << "\n\n";

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

void gsShapeOptProblem::runOptimization_aggregatedConstraints(index_t maxiter)
{
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
        *m_log << " Min d before updating: " << m_dJC->evalSubCon().minCoeff() << "\n\n";
        *m_log << " Constraint before updating: " << m_dJC->evalCon().minCoeff() << "\n\n";
        m_paramMethod->updateAndReset();
        m_paramMethod->computeMap();

        // FIXIT: log whether param was succesful
        *m_log << " Objective function after updating: " << evalObj() << "\n";
        *m_log << " Min d after updating: " << m_dJC->evalSubCon().minCoeff() << "\n\n";
        *m_log << " Constraint before updating: " << m_dJC->evalCon().minCoeff() << "\n\n";

        std::string nameAU = "cps_afterUpdate";
        m_log->saveVec(getFlat(),nameAU ,counter2);

        // gsMultiPatch<> dJ = m_dJC.getDetJ();
        // std::string namedj = "detJ";
        // m_log->plotMultiPatchOnGeometry(*m_mp,dJ,namedj,counter2);

        // namedj = "detJ_act";
        // m_log->plotActiveMultiPatchOnGeometry(*m_mp,dJ,m_dJC.eps(),namedj,counter2);

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
    // gsInfo << "obj: " << evalObj() << "\n";
    real_t obj = evalObj();
    // real_t winslow = m_paramMethod->evalObj();
    // real_t gradn = gradObj().norm();
    // real_t iters = (dynamic_cast< gsOptParamMethod *>(m_paramMethod))->iterations();

    gsVector<> v(1);
    v << obj;

    // std::string name = "geo";
    // m_log->plotInParaview(*m_mp,name,counter1);

    if (m_log->saveCps()){
        std::string name = "cps";
        m_log->saveVec(getFlat(),name,counter2,counter1++);
        m_log->logObj(v); // Recalculated, is there a better way?

        if (m_useOwnCons)
            *m_log << m_dJC->evalSubCon().minCoeff() <<  " " << m_dJC->evalCon().minCoeff() << "\n\n";
    }

    return true;
}

void gsShapeOptProblem::print()
{
  gsInfo << "m_numDesignVars  = " <<  m_numDesignVars << "\n";
  gsInfo << "m_numConstraints  = " <<  m_numConstraints << "\n";

  gsInfo << "m_eps = " << m_dJC->m_eps << "\n";

  gsInfo << "m_numConJacNonZero  = " <<  m_numConJacNonZero << "\n";
  gsInfo << "m_conJacRows.size  = " <<  m_conJacRows.size() << "\n";

  // for(std::vector<index_t>::iterator it = m_conJacRows.begin(); it != m_conJacCols.end(); ++it){
  //     gsInfo << " " << *it;
  // }

  gsInfo << "\nm_conJacCols.size  = " <<  m_conJacCols.size() << "\n";

  gsInfo << "m_conUpperBounds.size() = " << m_conUpperBounds.size() << "\n";
  gsInfo << "m_conLowerBounds.size() = " << m_conLowerBounds.size() << "\n";

  gsInfo << "\nm_conUpperBounds.max() = " << m_conUpperBounds.maxCoeff() << "\n";
  gsInfo << "m_conUpperBounds.min() = " << m_conUpperBounds.minCoeff() << "\n";

  gsInfo << "\nm_conLowerBounds.max() = " << m_conLowerBounds.maxCoeff() << "\n";
  gsInfo << "m_conLowerBounds.min() = " << m_conLowerBounds.minCoeff() << "\n";

  gsInfo << "m_desUpperBounds.size() = " << m_desUpperBounds.size() << "\n";
  gsInfo << "m_desLowerBounds.size() = " << m_desLowerBounds.size() << "\n";

  gsInfo << " \nMin d" << m_dJC->evalCon().minCoeff() << "\n";

  gsInfo << "m_curDesign.size() = " << m_curDesign.size() << "\n";

  gsMatrix<> disp(m_desUpperBounds.size(),3);
  disp << m_desLowerBounds,m_curDesign,m_desUpperBounds;
  // gsInfo << ".. design upper and lower bounds\n";
  // gsInfo << disp << "\n";
  //
  gsMatrix<> disp2(m_conUpperBounds.size(),2);
  disp2 << m_conLowerBounds,m_conUpperBounds;
  // gsInfo << ".. constraint upper and lower bounds\n";
  // gsInfo << disp2 << "\n";

}

void gsShapeOptProblem::setOptParamQuad(real_t quA, index_t quB)
{
    (dynamic_cast< gsAffineOptParamMethod *>(m_paramMethod.get()))->m_optParamMethod->setQuad(quA,quB);
}
