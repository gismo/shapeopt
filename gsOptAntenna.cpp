#include <gismo.h>
#include "gsOptAntenna.h"

// Implement
gsOptAntenna::gsOptAntenna(memory::shared_ptr<gsMultiPatch<>> mp, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, index_t param, real_t quA, index_t quB, bool useDetJCons, bool use_Lagrangian):
    gsShapeOptProblem(mp,slog,useDetJCons), m_stateEq(mp,numRefine)
{
    setupMappers();

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

    // Set the weight used in the objective function
    delta = memory::make_shared(new gsFunctionExpr<>("exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));
    ddeltadx = memory::make_shared(new gsFunctionExpr<>("-x/(0.1^2)*exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));
    ddeltady = memory::make_shared(new gsFunctionExpr<>("-y/(0.1^2)*exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));

}

gsOptAntenna::gsOptAntenna(memory::shared_ptr<gsMultiPatch<>> mp, index_t numRefine, memory::shared_ptr<gsShapeOptLog> slog, memory::shared_ptr<gsConstraint> constraint, index_t param, real_t quA, index_t quB, bool use_Lagrangian):
    gsShapeOptProblem(mp,slog,constraint), m_stateEq(mp,numRefine)
{

    if (use_Lagrangian)
        GISMO_ERROR("You cannot use Lagrangian for linearizations with aggregated constraints, hessian is not implemented yet");

    setupMappers();

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
    } else {
    // Setup optimization setupOptParameters
    // Calls the setupDesignBounds method
        GISMO_ERROR("Param value is not 0 or 1.. wrong input..\n");
    }

    // Setup optimization setupOptParameters
    // Calls the setupDesignBounds method
    setupOptParameters();

    // Set the weight used in the objective function
    delta = memory::make_shared(new gsFunctionExpr<>("exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));
    ddeltadx = memory::make_shared(new gsFunctionExpr<>("-x/(0.1^2)*exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));
    ddeltady = memory::make_shared(new gsFunctionExpr<>("-y/(0.1^2)*exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));

}

gsDofMapper gsOptAntenna::mapper_grad() const {
    if (m_mapper_grad_exist) // If they exist
        return m_mapper_grad;

    GISMO_ERROR("m_mapper_grad does not exist!");
}

real_t gsOptAntenna::evalObj() const {
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

    // variable deltaf = A.getCoeff(delta,G);
    variable df = A.getCoeff(*delta,G);

    // gsInfo<<"Plotting in Paraview...\n";
    // ev.options().setSwitch("plot.elements", true);
    // ev.writeParaview( deltaf    , G, "deltaf");
    variable u_r = A.getCoeff(u_real);
    variable u_i = A.getCoeff(u_imag);

    // return -2*ev.integral(df*meas(G));
    // return -2*ev.integral(df*u_r*meas(G));
    return -2*ev.integral(df*(u_r.sqNorm() + u_i.sqNorm())*meas(G));

}

gsVector<> gsOptAntenna::gradientObjWithoutAdjoint() const{
    gsMultiPatch<> u_real,u_imag;
    m_stateEq.solve(u_real,u_imag);

    gsMatrix<> dcdx = jacobDesignUpdate();
    gsMatrix<> dEdc = -2*(evaluateDerivativeTerm1(u_real,u_imag) + evaluateDerivativeTerm2(u_real,u_imag));

    gsMatrix<> out = dEdc.transpose()*dcdx;

    // gsInfo << "(" << out.rows() <<", " << out.cols() << ")\n" << std::flush;
    return out.transpose();
};

gsVector<> gsOptAntenna::gradAll() const{
    gsMultiPatch<> u_real,u_imag;
    m_stateEq.solve(u_real,u_imag);

    gsMatrix<> mat = m_stateEq.getDerivativeWithoutSolving(u_real,u_imag);

    gsVector<> dJdu = getObjDerivativeDu(u_real,u_imag);

    gsVector<> adjoint = m_stateEq.solveAdjoint(dJdu);

    // m_stateEq.printMatSize(mat,"mat");
    // m_stateEq.printMatSize(adjoint,"adjoint");
    gsVector<> term2 = adjoint.transpose()*mat;

    return -2*(evaluateDerivativeTerm1(u_real,u_imag) + term2);
}

// Adjoint method for sensitivities.
gsVector<> gsOptAntenna::gradObj() const{

    gsMatrix<> dcdx = jacobDesignUpdate();

    gsMatrix<> dEdc = mapGradient(m_mapper_grad, gradAll());

    gsMatrix<> out = dEdc.transpose()*dcdx;

    // gsInfo << "(" << out.rows() <<", " << out.cols() << ")\n" << std::flush;
    return out.transpose();
};

gsVector<> gsOptAntenna::evaluateDerivativeTerm1(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
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

    space u = A.getSpace(dbasis);

    variable u_r = A.getCoeff(u_real);
    variable u_i = A.getCoeff(u_imag);

    variable df = A.getCoeff(*delta,G);
    variable ddfdx = A.getCoeff(*ddeltadx,G);
    variable ddfdy = A.getCoeff(*ddeltady,G);

    A.initSystem();

    // Save the mapper used to calculate these guys..
    if (! m_mapper_grad_exist){
        m_mapper_grad = u.mapper();
        m_mapper_grad_exist = true;
    }


    gsFunctionExpr<> x("x",2);
    gsFunctionExpr<> y("y",2);
    variable fx = A.getCoeff(x);
    variable fy = A.getCoeff(y);
    auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
    auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
    auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
    auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

    auto g11 = j00*j00 + j10*j10;
    auto g12 = j00*j01 + j10*j11;
    auto g22 = j01*j01 + j11*j11;

    auto uxi = grad(u)*fjac(fx);
    auto ueta = grad(u)*fjac(fy);

    auto detJ = j00*j11 - j10*j01;
    auto detJinv = 1/detJ.val();

    auto signOfDetJ = detJinv*meas(G);

    auto d_detJ_dcx = signOfDetJ*(uxi*j11 - ueta*j10) ;
    auto d_detJ_dcy = signOfDetJ*(ueta*j00 - uxi*j01) ;

    auto term_1x = df.val()*d_detJ_dcx*(u_r.sqNorm() + u_i.sqNorm());
    auto term_2x = u*ddfdx.val()*(u_r.sqNorm() + u_i.sqNorm())*meas(G);
    // auto term_1x = df.val()*d_detJ_dcx*u_r;
    // auto term_2x = u*ddfdx.val()*u_r;

    A.assemble(term_1x + term_2x);

    gsMatrix<> xVec = A.rhs();

    A.initSystem();

    auto term_1y = df.val()*d_detJ_dcy*(u_r.sqNorm() + u_i.sqNorm());
    auto term_2y = u*ddfdy.val()*(u_r.sqNorm() + u_i.sqNorm())*meas(G);
    // auto term_1y = df.val()*d_detJ_dcy*u_r;
    // auto term_2y = u*ddfdy.val()*u_r;

    A.assemble(term_1y + term_2y);

    gsMatrix<> yVec = A.rhs();

    gsMatrix<> out(xVec.rows()*2,1);
    out << xVec, yVec;

    return out;
}

// Only used for gradientObjWithoutAdjoint
gsVector<> gsOptAntenna::evaluateDerivativeTerm2(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
    // Perhaps the solve is called unecessary in here: FIXIT
    gsMatrix<> dudc = m_stateEq.getDerivativeOfU();

    gsVector<> dJdu = getObjDerivativeDu(u_real,u_imag);

    // gsInfo << "\n\n-------RETURN SOME STUFF-------\n\n";
    // m_stateEq.printMatSize(dudc,"dudc");
    // m_stateEq.printMatSize(dJdu,"dJdu");
    return dudc*dJdu;
}

gsVector<> gsOptAntenna::getObjDerivativeDu(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
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

    variable df = A.getCoeff(*delta,G);

    A.initSystem();
    // A.assemble(df.val()*du*meas(G));
    A.assemble(df.val()*2*u_r.val()*du*meas(G));

    gsVector<> vec_re = A.rhs();

    A.initSystem();
    A.assemble(df.val()*2*u_i.val()*du*meas(G));
    // A.assemble(df.val()*2*u_i.val()*du*meas(G));

    gsVector<> vec_im = A.rhs();

    gsVector<> out;
    out.setZero(A.rhs().rows()*2);
    out << vec_re,vec_im;

    return out;
}

void gsOptAntenna::setupMappers()
{
    // Get mappers from multibasis with interfaces glued
    gsMultiBasis<> geoBasis(*m_mp);
    geoBasis.getMapper(iFace::glue,m_mappers[0],false); // False means that we do not finalize
    geoBasis.getMapper(iFace::glue,m_mappers[1],false); // False means that we do not finalize

    // Fix y coordinate for all boundaries
    for (index_t i = 0; i < m_mp->nBoundary(); i ++){
        patchSide ps = m_mp->boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
        // gsInfo << ps << " is fixed in y direction \n\n";
        // gsInfo << "boundaryDofs on patch " << ps.patch << "\n"<< boundaryDofs << "\n\n";
        m_mappers[1].markBoundary(ps.patch,boundaryDofs);

        if (m_mp->nBoxes() > 3){
            if (ps.patch == 4){
                m_mappers[0].markBoundary(ps.patch,boundaryDofs);
            }
        } else {
            // If there is less than 3 patches then all boundaries is fixed
            m_mappers[0].markBoundary(ps.patch,boundaryDofs);
        }
    }

    // Fix boundary of the fixedPatch
    gsMultiPatch<> fixedGeom(m_mp->patch(fixedPatch));

    for(index_t i = 0; i < fixedGeom.nBoundary(); i++){
        patchSide ps = fixedGeom.boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(fixedPatch).boundary(ps);
        m_mappers[0].markBoundary(fixedPatch,boundaryDofs); // Mark xcoord of boundaries of fixed patch
        m_mappers[1].markBoundary(fixedPatch,boundaryDofs); // Mark ycoord of boundaries of
    }

    // Finalize mappers
    m_mappers[0].finalize();
    m_mappers[1].finalize();

    gsInfo << "n controlpoints: " << m_mappers[0].mapSize() << "\n";

    // Tag fixed bnds
    for(index_t i = 0; i < fixedGeom.nBoundary(); i++){
        patchSide ps = fixedGeom.boundaries()[i];
        gsVector<unsigned> boundaryDofs = m_mp->basis(fixedPatch).boundary(ps);
        for (index_t j = 0; j < boundaryDofs.size(); j ++){
            m_mappers[0].markTagged(boundaryDofs[j],fixedPatch); // Mark xcoord of boundaries of fixed patch
            m_mappers[1].markTagged(boundaryDofs[j],fixedPatch); // Mark ycoord of boundaries of

        }
    }

    // count number of free dofs
    n_free = 0;
    for(index_t i = 0; i < m_mp->targetDim(); i++){
        n_free += m_mappers[i].freeSize();
    }
}

void gsOptAntenna::setupDesignBounds()
{

    // Set design bounds
    real_t lb_x = -d_bw*m_stateEq.pde_L_f/2;
    real_t ub_x = d_bw*m_stateEq.pde_L_f/2;
    real_t lb_y = d_g/2*m_stateEq.pde_L_f;
    real_t ub_y = (d_g/2 + d_bh)*m_stateEq.pde_L_f;

    gsInfo << "width = " << (ub_x - lb_x)/m_stateEq.pde_r_t << " \n";
    gsInfo << "height = " << (ub_y - lb_y)/m_stateEq.pde_r_t << " \n";
    gsInfo << "gap = " << (lb_y)/m_stateEq.pde_r_t << " \n";

    m_desLowerBounds.setZero(n_tagged);
    // X coordinates
    gsVector<> xLowerBounds;
    xLowerBounds.setConstant(m_mappers[0].taggedSize(),lb_x);
    m_desLowerBounds.segment(0,m_mappers[0].taggedSize()) = xLowerBounds;

    // Y coordinates
    gsVector<> yLowerBounds;
    yLowerBounds.setConstant(m_mappers[1].taggedSize(),lb_y);
    m_desLowerBounds.segment(m_mappers[0].taggedSize(),m_mappers[1].taggedSize()) = yLowerBounds;

    m_desUpperBounds.setZero(n_tagged);
    // X coordinates
    gsVector<> xUpperBounds;
    xUpperBounds.setConstant(m_mappers[0].taggedSize(),ub_x);
    m_desUpperBounds.segment(0,m_mappers[0].taggedSize()) = -xLowerBounds;

    // Y coordinates
    gsVector<> yUpperBounds;
    yUpperBounds.setConstant(m_mappers[1].taggedSize(),ub_y);
    m_desUpperBounds.segment(m_mappers[0].taggedSize(),m_mappers[1].taggedSize()) = yUpperBounds;
}
