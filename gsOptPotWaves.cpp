#include <gismo.h>
#include "gsOptPotWaves.h"
#include "gsMyExpressions.h"

// Implement
gsOptPotWaves::gsOptPotWaves(gsMultiPatch<>::Ptr mp, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, index_t param, real_t quA, index_t quB, bool useDetJCons, bool use_Lagrangian):
    gsShapeOptProblem(mp,slog,useDetJCons), m_stateEq(mp, numRefine)
{
    setupMappers();

    constructor(param, quA, quB, use_Lagrangian);

}

void gsOptPotWaves::constructor(index_t param, real_t quA, index_t quB, bool use_Lagrangian)
{
    *m_log << "quA for parametrization: " << quA << "\n";
    *m_log << "quB for parametrization: " << quB << "\n\n";

    // Allocate parametrization method
    if (param == 0){
        m_paramMethod = memory::make_shared(new gsSpringMethod(m_mp,m_mappers));
        m_paramMethod->computeMap();
    } else if (param == 1) {
        gsModLiao::Ptr opt_param = memory::make_shared(new gsModLiao(m_mp,m_mappers,true));
        opt_param->setQuad(quA,quB);
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 2) {
        gsWinslow::Ptr opt_param = memory::make_shared(new gsWinslow(m_mp,m_mappers,true));
        opt_param->setQuad(quA,quB);// FIXIT: this and the next two statements can be moved out of if state ment if param =! 0 if ()...
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 3) {
        gsLiao::Ptr opt_param = memory::make_shared(new gsLiao(m_mp,m_mappers,true));
        opt_param->setQuad(quA,quB);
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 4) {
        gsHarmonic::Ptr opt_param = memory::make_shared(new gsHarmonic(m_mp,m_mappers,true));
        opt_param->setQuad(quA,quB);
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else if (param == 5) {
        m_paramMethod = memory::make_shared(new gsWinslowWithDeriv(m_mp,m_mappers,false,false,true,0));
        // m_paramMethod->setQuad(quA,quB); FIXIT: implement, maybe with dynamic cast
    } else if (param == 6) {
        gsWinslow::Ptr opt_param = memory::make_shared(new gsWinslow(m_mp,m_mappers,false,false,true,0));
        // opt_param->setQuad(quA,quB);// FIXIT: this and the next two statements can be moved out of if state ment if param =! 0 if ()...
        opt_param->setQuad(4,4);
        *m_log << "quA, quB = " << opt_param->m_quA << ", " << opt_param->m_quB << "\n";
        m_paramMethod = memory::make_shared(new gsAffineOptParamMethod(opt_param, use_Lagrangian));
        m_paramMethod->computeMap();
    } else {
        GISMO_ERROR("Param value is not 0 to 6... wrong input..\n");
    }

    // Copy some data from m_paramMethod for ease of use
    n_free = m_paramMethod->n_free;
    n_flat = m_paramMethod->n_flat;
    n_tagged = m_paramMethod->n_tagged;
    n_cps = m_paramMethod->n_cps;

    // Setup optimization setupOptParameters
    // Calls the setupDesignBounds method
    setupOptParameters();
}

gsDofMapper gsOptPotWaves::mapper_grad() const {
    if (m_mapper_grad_exist) // If they exist
        return m_mapper_grad;

    GISMO_ERROR("m_mapper_grad does not exist!");
}

real_t gsOptPotWaves::evalObj() const {
    gsMultiPatch<> u_real,u_imag;
    m_stateEq.solve(u_real,u_imag);

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    ev.options().setInt("quB",m_stateEq.quA()); // FIXIT lower no points
    ev.options().setReal("quA",m_stateEq.quB()); // FIXIT lower no points

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);
    gsMultiBasis<> dbasis(m_stateEq.dbasis);
    A.setIntegrationElements(dbasis);

    // gsInfo<<"Plotting in Paraview...\n";
    // ev.options().setSwitch("plot.elements", true);
    // ev.writeParaview( deltaf    , G, "deltaf");
    variable u_r = A.getCoeff(u_real);
    variable u_i = A.getCoeff(u_imag);
    //

    // FIXIT: shouldn't I multiply with 0.5 for this to work?
    return ev.integral(u_r.sqNorm());

}

gsVector<> gsOptPotWaves::gradAll() const{
    gsMultiPatch<> u_real,u_imag;
    m_stateEq.solve(u_real,u_imag);

    gsMatrix<> mat = m_stateEq.getDerivativeWithoutSolving(u_real,u_imag);

    gsVector<> dJdu = getObjDerivativeDu(u_real,u_imag);

    gsVector<> adjoint = m_stateEq.solveAdjoint(dJdu);

    // m_stateEq.printMatSize(mat,"mat");
    // m_stateEq.printMatSize(adjoint,"adjoint");
    gsVector<> term2 = adjoint.transpose()*mat;

    return (getObjDerivativeDc(u_real,u_imag) + term2);
}

// Adjoint method for sensitivities.
gsVector<> gsOptPotWaves::gradObj() const{
    // get gradObj from gradAll!
    
    gsMatrix<> dcdx = jacobDesignUpdate();

    gsMatrix<> dEdc = mapGradient(m_mapper_grad, gradAll());

    gsMatrix<> out = dEdc.transpose()*dcdx;

    // gsInfo << "(" << out.rows() <<", " << out.cols() << ")\n" << std::flush;
    return out.transpose();
   
};

void gsOptPotWaves::setupMapperGrad()
{    
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(m_stateEq.dbasis);

    typedef gsExprAssembler<>::space       space;

    space u = A.getSpace(dbasis);

    A.initSystem();

    m_mapper_grad = u.mapper();
}

void gsOptPotWaves::setupMappers()
{
    // Get mappers from multibasis with interfaces glued
    gsMultiBasis<> geoBasis(*m_mp);

    for( index_t d = 0; d < m_dim; d++)
        geoBasis.getMapper(iFace::glue,m_mappers[d],false); // False means that we do not finalize

    // Coordinates of all boundaries
    for (index_t i = 0; i < m_mp->nBoundary(); i++){
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
        // gsInfo << ps << " is fixed in y direction \n\n";
        // gsInfo << "boundaryDofs on patch " << ps.patch << "\n"<< boundaryDofs << "\n\n";
        for( index_t d = 0; d < m_dim; d++)
            m_mappers[d].markBoundary(ps.patch,boundaryDofs);
    }


    gsInfo << "n controlpoints: " << m_mappers[0].mapSize() << "\n";
    n_controlpoints = m_mappers[0].mapSize() ;

    // FIXIT : Fix all non domain patches
    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        if (!isPatchInDomain(p)) 
        {
            eliminatePatch(p);
        }
    }

    // Finalize mappers
    for( index_t d = 0; d < m_dim; d++)
        m_mappers[d].finalize();


    // Tag Gamma_S
    for(index_t i = 0; i < m_mp->nBoundary(); i++){
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

        if (m_stateEq.isBndGamma_s(i))
        {
            for (index_t j = 0; j < boundaryDofs.size(); j ++)
            {
                for( index_t d = 0; d < m_dim; d++)
                    m_mappers[d].markTagged(boundaryDofs[j],ps.patch); 
            }
        }
    }

    // count number of free dofs
    n_free = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_free += m_mappers[i].freeSize();
    }

    // Setup mapper_grad
    setupMapperGrad();
}

void gsOptPotWaves::eliminatePatch(index_t p)
{
    for (index_t i = 0; i < m_mp->patch(p).coefsSize(); i++)
    {
        for (index_t d = 0; d < m_dim; d++)
        {
            m_mappers[d].eliminateDof(i, p);
        }
    }
}

//FIXIT: implement
void gsOptPotWaves::setupDesignBounds()
{
    m_desLowerBounds.setZero(n_tagged);
    m_desUpperBounds.setZero(n_tagged);
}

// Derivative of obj wrt c (not inlcuding du/dc terms)
gsVector<> gsOptPotWaves::getObjDerivativeDc(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(m_stateEq.dbasis);
    gsExprEvaluator<> ev(A);

    A.options().setReal("quA",m_stateEq.quA());
    A.options().setInt("quB",m_stateEq.quB());

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    space u = A.getSpace(dbasis,m_dim);

    variable u_r = A.getCoeff(u_real);
    variable u_i = A.getCoeff(u_imag);

    A.initSystem();

    gsMatrix<> out;
    out.setZero(A.rhs().rows(),1);

    return out;
}

gsVector<> gsOptPotWaves::getObjDerivativeDu(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(m_stateEq.dbasis);
    gsExprEvaluator<> ev(A);

    A.options().setReal("quA",m_stateEq.quA());
    A.options().setInt("quB",m_stateEq.quB());

    // opts.setReal("quA",2);
    // gsOptionList opts = A.options();
    // opts.setInt("quB",1);
    // gsInfo << opts << "\n";
    // A.setOptions(opts);

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    space du = A.getSpace(m_stateEq.dbasis);
    du.setInterfaceCont(0);

    variable u_r = A.getCoeff(u_real);
    variable u_i = A.getCoeff(u_imag);

    A.initSystem();
    // A.assemble(df.val()*du*meas(G));
    A.assemble(2*du*u_r);

    gsVector<> vec_re = A.rhs();

    A.initSystem();
    //A.assemble(df.val()*2*u_i.val()*du*meas(G));
    // A.assemble(df.val()*2*u_i.val()*du*meas(G));

    gsVector<> vec_im = A.rhs();

    gsVector<> out;
    out.setZero(A.rhs().rows()*2);
    out << vec_re,vec_im;

    return out;
}

bool gsOptPotWaves::isFlatInPML(index_t i)
{
    GISMO_NO_IMPLEMENTATION;
}

bool gsOptPotWaves::isPatchInDomain(index_t p)
{
    return m_stateEq.isPatchInDomain(p);
}

bool gsOptPotWaves::isCpsInDomain(index_t i, index_t p, index_t dim)
{
    // If (i,p) is not coupled, just return whether p is a PML patch
    if (!m_mappers[dim].is_coupled(i,p))
        return m_stateEq.isPatchInDomain(p);

    // Else we need to find the coupled patches
    index_t ii = m_mappers[dim].index(i,p);
    std::vector< std::pair< index_t, index_t > > result;

    m_mappers[dim].preImage(ii, result);

    for(std::vector<std::pair<index_t,index_t> >::const_iterator it = result.begin(); it != result.end(); ++it)
   {
       // Get local patch
       index_t patch = it->first;

       if ( ! isPatchInDomain(patch))
           return false;

   }

}

