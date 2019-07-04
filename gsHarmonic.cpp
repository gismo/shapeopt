#include <gismo.h>
#include "gsHarmonic.h"
using namespace gismo;

real_t gsHarmonic::evalObj() const {
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    gsExprEvaluator<> ev(A);

    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);

    A.setIntegrationElements(dbasis);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);
    gsFunctionExpr<> f1("x", "-y",2);
    variable ffun = ev.getVariable(f1);
    auto D = fjac(ffun);

    gsMultiPatch<> x_mp = m_mp->coord(0);
    gsMultiPatch<> y_mp = m_mp->coord(1);

    variable x = A.getCoeff(x_mp);
    variable y = A.getCoeff(y_mp);

    // auto Lx = fform(G).adj() % hess(x); // Why does this not work?
    auto Lx = (jac(G).tr()*jac(G)).adj() % hess(x);
    // auto Ly = fform(G).adj() % hess(y);
    auto Ly = (jac(G).tr()*jac(G)).adj() % hess(y);

    real_t out = 0;

    out += ev.integral(lambda_1 * hess(G).sqNorm());
    out += ev.integral(lambda_2 * jac(G).sqNorm());
    // out += ev.integral( (D*hess(x)*D*fform(G)).trace().sqNorm().val() );
    // out += ev.integral( (D*hess(y)*D*fform(G)).trace().sqNorm().val() );
    out += ev.integral(Lx.sqr());
    out += ev.integral(Ly.sqr());
    return out;
}

// gsVector<> gsHarmonic::gradObj() const{
//     gsVector<> u = getFree();
//     const index_t n = u.rows();
//     //GISMO_ASSERT((index_t)m_numDesignVars == n*m, "Wrong design.");
//
//     gsVector<> result(n);
//
//     gsMatrix<real_t> uu = u;//copy
//     gsAsVector<real_t> tmp(uu.data(), n);
//     gsAsConstVector<real_t> ctmp(uu.data(), n);
//     index_t c = 0;
//
//     // for all partial derivatives (column-wise)
//     for ( index_t i = 0; i!=n; i++ )
//     {
//         // to do: add m_desLowerBounds m_desUpperBounds check
//         tmp[i]  += real_t(0.00001);
//         updateFree(tmp);
//         const real_t e1 = this->evalObj();
//         tmp[i]   = u[i] + real_t(0.00002);
//         updateFree(tmp);
//         const real_t e3 = this->evalObj();
//         tmp[i]   = u[i] - real_t(0.00001);
//         updateFree(tmp);
//         const real_t e2 = this->evalObj();
//         tmp[i]   = u[i] - real_t(0.00002);
//         updateFree(tmp);
//         const real_t e4 = this->evalObj();
//         tmp[i]   = u[i];
//         updateFree(tmp);
//         result[c++]= ( 8 * (e1 - e2) + e4 - e3 ) / real_t(0.00012);
//     }
//     return result;
// }

// FIXIT: make dimension independent
gsVector<> gsHarmonic::gradObj() const{
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);

    gsOptionList opts = A.options();
    opts.setInt("quB",m_quB);
    opts.setReal("quA",m_quA);
    A.setOptions(opts);

    A.setIntegrationElements(dbasis);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    space u = A.getSpace(dbasis);

    gsMultiPatch<> x_mp = m_mp->coord(0);
    gsMultiPatch<> y_mp = m_mp->coord(1);

    variable x = A.getCoeff(x_mp);
    variable y = A.getCoeff(y_mp);

    auto Lx = (jac(G).tr()*jac(G)).adj() % hess(x);
    // auto Lx = fform(G).adj() % hess(x);
    // auto Ly = fform(G).adj() % hess(y);
    auto Ly = (jac(G).tr()*jac(G)).adj() % hess(y);

    auto dLxdcx = 2*(jac(G).tr()*dJacdc(u,0)).adj()%hess(x) + hess(u) %(jac(G).tr()*jac(G)).adj();
    auto dLxdcy = 2*(jac(G).tr()*dJacdc(u,1)).adj()%hess(x);

    auto dLydcx = 2*(jac(G).tr()*dJacdc(u,0)).adj()%hess(y);
    auto dLydcy = 2*(jac(G).tr()*dJacdc(u,1)).adj()%hess(y) + hess(u) %(jac(G).tr()*jac(G)).adj();

    A.initSystem();

    // FIXIT: make dimension independent
    // x coordinate
    A.assemble(2*lambda_1*hess(u)%hess(x));
    A.assemble(2*lambda_2*dJacdc(u,0)%jac(G));
    A.assemble(2*Lx.val()*dLxdcx);
    A.assemble(2*Ly.val()*dLydcx);
    gsVector<> xVec = A.rhs();

    A.initSystem();

    // y coordinate

    A.assemble(2*lambda_1*hess(u)%hess(y));
    A.assemble(2*lambda_2*dJacdc(u,1)%jac(G));
    A.assemble(2*Lx.val()*dLxdcy);
    A.assemble(2*Ly.val()*dLydcy);
    gsVector<> yVec = A.rhs();


    gsVector<> all(m_mp->targetDim()*dbasis.size());
    all << xVec, yVec;
    // all += gradObjHelper();
    return mapMatrix(u.mapper(),all);
}

gsMatrix<> gsHarmonic::hessObj(gsMatrix<> &hessObjTagged) const{
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);

    gsOptionList opts = A.options();
    opts.setInt("quB",m_quB);
    opts.setReal("quA",m_quA);
    A.setOptions(opts);

    A.setIntegrationElements(dbasis);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    space u = A.getSpace(dbasis);

    A.initSystem();

    gsFunctionExpr<> x("x",2);
    gsFunctionExpr<> y("y",2);
    variable ffx = A.getCoeff(x);
    variable ffy = A.getCoeff(y);
    auto j00 = fjac(ffx).tr()*jac(G)*fjac(ffx);
    auto j10 = fjac(ffy).tr()*jac(G)*fjac(ffx);
    auto j01 = fjac(ffx).tr()*jac(G)*fjac(ffy);
    auto j11 = fjac(ffy).tr()*jac(G)*fjac(ffy);

    auto g11 = j00*j00 + j10*j10;
    auto g12 = j00*j01 + j10*j11;
    auto g22 = j01*j01 + j11*j11;

    auto uxi = grad(u)*fjac(ffx);
    auto ueta = grad(u)*fjac(ffy);

    auto detJ = j00*j11 - j01*j10;
    auto detJinv = detJ.inv();
    auto detJinv2 = detJinv*detJinv;
        auto detJinv3 = detJinv2*detJinv;

      auto d_detJ_dcx = uxi*j11 - ueta*j10 ;
      auto d_g11pg22_dcx = 2*(uxi*j00 + ueta*j01);

      auto fx = -d_detJ_dcx*(g11 + g22);
      auto gx = d_g11pg22_dcx;

      auto I = (g11+g22)*detJinv;

      auto d_I_dcx = fx*detJinv2 + gx*detJinv;

      auto d_detJinv2_dcx = -2 * d_detJ_dcx * detJinv3;
      auto d_detJinv_dcx = - d_detJ_dcx * detJinv2;

      auto d_fx_dcx_detJinv2 = -d_detJ_dcx*detJinv2*d_g11pg22_dcx.tr();
      auto d_gx_dcx_detJinv = 2*(uxi*detJinv*uxi.tr() + ueta*detJinv*ueta.tr());

      auto d2_I_dcx2 = fx * d_detJinv2_dcx.tr() + d_fx_dcx_detJinv2 + gx*d_detJinv_dcx.tr() + d_gx_dcx_detJinv;

      A.assemble(-d2_I_dcx2);

      gsMatrix<> xxMat = A.matrix();

      A.initSystem();

      auto d_detJ_dcy = ueta*j00 - uxi*j01 ;
      auto d_g11pg22_dcy = 2*(uxi*j10 + ueta*j11);

      auto fy = -d_detJ_dcy*(g11 + g22);
      auto gy = d_g11pg22_dcy;


      auto d_I_dcy = fy*detJinv2 + gy*detJinv;

      auto d_detJinv2_dcy = -2 * d_detJ_dcy * detJinv3;
      auto d_detJinv_dcy = - d_detJ_dcy * detJinv2;

      auto d_fx_dcy_detJinv2 = ueta*(g11 + g22)*detJinv2*uxi.tr() - uxi*(g11 + g22)*detJinv2*ueta.tr() - d_detJ_dcx*detJinv2*d_g11pg22_dcy.tr();

      auto d2_I_dcxcy = fx * d_detJinv2_dcy.tr() + d_fx_dcy_detJinv2 + gx*d_detJinv_dcy.tr();

      A.assemble(-d2_I_dcxcy);
      gsMatrix<> xyMat = A.matrix().transpose();

      A.initSystem();

      auto d_fy_dcy_detJinv2 = -d_detJ_dcy*detJinv2*d_g11pg22_dcy.tr();
      auto d_gy_dcy_detJinv = 2*(uxi*detJinv*uxi.tr() + ueta*detJinv*ueta.tr());

      auto d2_I_dcy2 = fy*d_detJinv2_dcy.tr() + d_fy_dcy_detJinv2 + gy*d_detJinv_dcy.tr() + d_gy_dcy_detJinv;

      A.assemble(-d2_I_dcy2);


      gsMatrix<> yyMat = A.matrix();

        gsMatrix<> all(xxMat.rows() + yyMat.rows(),xxMat.cols() + yyMat.cols());

        all << xxMat, xyMat.transpose(), xyMat, yyMat;

        // Map twice, is there a better way?
        gsMatrix<> all2 = mapMatrix(u.mapper(),all);

        // Map tagged part
        hessObjTagged = mapMatrixToTagged(u.mapper(),all2);

        return mapMatrix(u.mapper(),all2);

    }
