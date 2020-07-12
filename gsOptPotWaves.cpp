#include <gismo.h>
#include "gsOptPotWaves.h"
#include "gsMyExpressions.h"
#include "gsStateEquationPotWaves.h"

// Implement
gsOptPotWaves::gsOptPotWaves(gsMultiPatch<>::Ptr mp, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, index_t param, real_t quA, index_t quB, bool useDetJCons, bool use_Lagrangian):
    gsShapeOptProblem(mp,slog,useDetJCons), m_stateEq(mp, numRefine)
{
    setupMappers();

    constructor(param, quA, quB, use_Lagrangian);

    getTestMatrix();

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
    //gsMultiPatch<> u_real,u_imag;
    //gsVector<> uvec_real, uvec_imag;
    //m_stateEq.solve(u_real,u_imag, uvec_real, uvec_imag);

    //getObjVec();
    //gsVector<> out =  m_objVec_re.transpose()*uvec_real + m_objVec_im.transpose()*uvec_imag;
    //return out[0];
    
    gsVector<> flat = m_paramMethod->getFlat();

    return 0.5*flat.transpose()*m_testMatrix*flat;

    /*
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
    gsDebugVar(out);
    gsDebugVar(ev.integral(u_r));
    return ev.integral(u_r);
    */
    

}

gsVector<> gsOptPotWaves::gradAll() const{
    //gsMultiPatch<> u_real,u_imag;
    //gsVector<> uvec_real, uvec_imag;
    //m_stateEq.solve(u_real,u_imag, uvec_real, uvec_imag);

    //gsMatrix<> mat = m_stateEq.getDerivativeWithoutSolving(u_real,u_imag);

    //gsVector<> dJdu = getObjDerivativeDu(u_real,u_imag);

    //gsVector<> adjoint = m_stateEq.solveAdjoint(dJdu);

    //// m_stateEq.printMatSize(mat,"mat");
    //// m_stateEq.printMatSize(adjoint,"adjoint");
    //gsVector<> term2 = adjoint.transpose()*mat;

    //return (getObjDerivativeDc(uvec_real,uvec_imag) + term2);
    //
    gsVector<> flat = m_paramMethod->getFlat();

    return m_testMatrix*flat;
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
    m_mapper_grad_exist = true;
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
                {
                    index_t i = boundaryDofs[j];

                    // Skip if cps(i,p,d) == 0
                    if (m_mp->patch(ps.patch).coef(i,d) == 0) continue;

                    m_mappers[d].markTagged(i,ps.patch); 
                }
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

//FIXIT: implemen
void gsOptPotWaves::setupDesignBounds()
{
    m_desLowerBounds.setConstant(n_tagged,-1e3);
    m_desUpperBounds.setConstant(n_tagged,1e3);

    gsVector<> tagged = m_paramMethod->getTagged();
    for (index_t i = 0; i < n_tagged; i++)
    {
        if( i < m_paramMethod->m_shift_tagged[1])     // Direction x
            continue;
        else if( i < m_paramMethod->m_shift_tagged[2])// Direction y
            m_desLowerBounds[i] = 0;
        else                            // Direction z
            m_desUpperBounds[i] = 0;
    }

    //gsMatrix<> disp(n_tagged,3);
    //disp << m_desLowerBounds, tagged, m_desUpperBounds;
    //gsInfo << "\n" << disp << "\n";
}

// Derivative of obj wrt c (not inlcuding du/dc terms)
gsVector<> gsOptPotWaves::getObjDerivativeDc(gsVector<> &u_real, gsVector<> &u_imag) const{

    getGradObjMat();

    return m_gradObjMat_re*u_real + m_gradObjMat_im*u_imag;

    /*
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
    */
}

gsVector<> gsOptPotWaves::getObjDerivativeDu(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const
{
    /*
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
    A.assemble(du);

    gsVector<> vec_re = A.rhs();

    A.initSystem();
    //A.assemble(df.val()*2*u_i.val()*du*meas(G));
    // A.assemble(df.val()*2*u_i.val()*du*meas(G));

    gsVector<> vec_im = A.rhs();
    */
    getObjVec();

    gsVector<> out;
    out.setZero(m_objVec_re.rows()*2);
    out << m_objVec_re,m_objVec_im;
    

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

void gsOptPotWaves::getObjVec() const {

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

    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);


    /*gsFunctionExpr<> zerfun("0.0",3);
    variable zer = ev.getVariable(zerfun);
    auto zero_matrix = zer.val()*u*u.tr();  // Perhaps I need to multiply with u*u.tr();
    auto zero_vec = zer.val()*u*nv(G).norm();     // Perhaps I need to multiply with u;
    */

    gsFunctionExpr<>::uPtr expKzpiKx_re, expKzpiKx_im;
    m_stateEq.getObjFunctions( expKzpiKx_re, expKzpiKx_im);

    variable exp_re = A.getCoeff(*expKzpiKx_re,G);
    variable exp_im = A.getCoeff(*expKzpiKx_im,G);

    gsFunctionExpr<>::uPtr grad_expKzpiKx_re, grad_expKzpiKx_im, hess_asVec_expKzpiKx_re, hess_asVec_expKzpiKx_im;
    m_stateEq.getObjFunctions(grad_expKzpiKx_re, grad_expKzpiKx_im, hess_asVec_expKzpiKx_re, hess_asVec_expKzpiKx_im);

    variable gexp_re = A.getCoeff(*grad_expKzpiKx_re,G);
    variable gexp_im = A.getCoeff(*grad_expKzpiKx_im,G);

    A.initSystem();
    A.assemble(exp_im.val()*igrad(u,G)*nv(G));//*nv(G).norm());
    A.assemble( - u*(gexp_im.tr()*nv(G)));//*nv(G).norm());
    m_objVec_re =  - 2 * A.rhs();

    A.initSystem();
    A.assemble(exp_re.val()*igrad(u,G)*nv(G));//*nv(G).norm());
    A.assemble( - u*(gexp_re.tr()*nv(G)));//*nv(G).norm());
    m_objVec_im =  - 2 * A.rhs();

}

void gsOptPotWaves::getGradObjMat() const
{
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
    gsMultiBasis<> gbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);

    space v = A.getTestSpace(u,gbasis,m_dim);

    /*gsFunctionExpr<> zerfun("0.0",3);
    variable zer = ev.getVariable(zerfun);
    auto zero_matrix = zer.val()*u*u.tr();  // Perhaps I need to multiply with u*u.tr();
    auto zero_vec = zer.val()*u*nv(G).norm();     // Perhaps I need to multiply with u;
    */

    gsFunctionExpr<>::uPtr expKzpiKx_re, expKzpiKx_im;
    m_stateEq.getObjFunctions( expKzpiKx_re, expKzpiKx_im);

    variable exp_re = A.getCoeff(*expKzpiKx_re,G);
    variable exp_im = A.getCoeff(*expKzpiKx_im,G);

    gsFunctionExpr<>::uPtr grad_expKzpiKx_re, grad_expKzpiKx_im, hess_asVec_expKzpiKx_re, hess_asVec_expKzpiKx_im;
    m_stateEq.getObjFunctions(grad_expKzpiKx_re, grad_expKzpiKx_im, hess_asVec_expKzpiKx_re, hess_asVec_expKzpiKx_im);

    variable gexp_re = A.getCoeff(*grad_expKzpiKx_re,G);
    variable gexp_im = A.getCoeff(*grad_expKzpiKx_im,G);

    variable h_asVec_exp_re = A.getCoeff(*hess_asVec_expKzpiKx_re,G);
    variable h_asVec_exp_im = A.getCoeff(*hess_asVec_expKzpiKx_im,G); 

    //gsDebugVar(*grad_expKzpiKx_im);
    //gsDebugVar(*hess_asVec_expKzpiKx_im);

    auto hess_re = reshape(h_asVec_exp_re,m_dim,m_dim);
    auto hess_im = reshape(h_asVec_exp_im,m_dim,m_dim);

    // Obs I had to multiply jac(G).inv() to the other term in 'collapse' to get the code to run
    auto dJinvTdc = - matrix_by_space_tr(jac(G).inv(),jac(v));//*jac(G).inv();

    // FIXIT: Use unit normal vec here?
    A.initSystem();

    //A.assemble(- (nvDeriv(v,G)*nv(G)/nv(G).norm())*u.tr()*(gexp_im.tr()*nv(G)).val());
    A.assemble(- (nvDeriv(v,G)*gexp_im)*u.tr());//*nv(G).norm());
    A.assemble(- (v*hess_im*nv(G))*u.tr());//*nv(G).norm());

    A.assemble((v*gexp_im)*(igrad(u,G)*nv(G)).tr());//*nv(G).norm());
    A.assemble(exp_im.val()*my_collapse(nv(G).tr(),dJinvTdc)*jac(G).inv().tr()*grad(u).tr());//*nv(G).norm());
    A.assemble(exp_im.val()*nvDeriv(v,G)*igrad(u,G).tr());//)*nv(G).norm() );
    //A.assemble((nvDeriv(v,G)*nv(G)/nv(G).norm())*exp_im.val()*(igrad(u,G)*nv(G)).tr() );

    m_gradObjMat_re =  - 2 * A.matrix();

    A.initSystem();
    //A.assemble(- (nvDeriv(v,G)*nv(G)/nv(G).norm())*u.tr()*(gexp_re.tr()*nv(G)).val());
    A.assemble(- (nvDeriv(v,G)*gexp_re)*u.tr());//*nv(G).norm());
    A.assemble(- (v*hess_re*nv(G))*u.tr());//*nv(G).norm());

    A.assemble((v*gexp_re)*(igrad(u,G)*nv(G)).tr());//*nv(G).norm());
    A.assemble(exp_re.val()*my_collapse(nv(G).tr(),dJinvTdc)*jac(G).inv().tr()*grad(u).tr());//*nv(G).norm());
    A.assemble(exp_re.val()*nvDeriv(v,G)*igrad(u,G).tr());//*nv(G).norm() );
    //A.assemble((nvDeriv(v,G)*nv(G)/nv(G).norm())*exp_re.val()*(igrad(u,G)*nv(G)).tr() );

    m_gradObjMat_im =  - 2 * A.matrix();

}

void gsOptPotWaves::getTestMatrix()
{
    m_testMatrix.setZero(n_flat,n_flat);
    index_t flat = 0;

    for(index_t d = 0; d < m_mp->targetDim(); ++d)
    {
        for (index_t p = 0; p < m_mp->nBoxes(); p++)
        {
            for (index_t i = 0; i < m_mp->patch(p).coefsSize(); i++)
            {

                if (m_mappers[d].is_tagged(i,p))
                {
                    m_testMatrix(flat,flat) = 1;
                }
                
                flat++;
            }
        }
    }

}
    
