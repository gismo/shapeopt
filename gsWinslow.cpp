#include <gismo.h>
#include "gsWinslow.h"
using namespace gismo;

real_t gsWinslow::evalObj() const {
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);
    auto detJinv = jac(G).inv().det(); // The inverse of det J

    return ev.integral((jac(G).tr()*jac(G)).trace()*detJinv);
}

gsVector<> gsWinslow::gradObj() const{
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

    space u = A.getSpace(dbasis,m_mp->geoDim()); // The gradient is a vector with targetDim values

    A.initSystem();

    auto detJinv = jac(G).inv().det(); // The inverse of abs(det J)
    auto detJinv2 = detJinv*detJinv;

    A.assemble(2*matrix_by_space(jac(G).tr(),jac(u)).trace()*detJinv
        - matrix_by_space(jac(G).inv(),jac(u)).trace()*
        (jac(G).tr()*jac(G)).trace().val()*detJinv);

    gsVector<> all = A.rhs();
    return mapMatrix(u.mapper(),all);

}

// gsMatrix<> gsWinslow::hessObj(gsMatrix<> &hessObjTagged) const{
//     gsExprAssembler<> A(1,1);
//     gsMultiBasis<> dbasis(*m_mp);
//
//     gsOptionList opts = A.options();
//     opts.setInt("quB",m_quB);
//     opts.setReal("quA",m_quA);
//     A.setOptions(opts);
//
//     A.setIntegrationElements(dbasis);
//
//     typedef gsExprAssembler<>::geometryMap geometryMap;
//     typedef gsExprAssembler<>::variable    variable;
//     typedef gsExprAssembler<>::space       space;
//     typedef gsExprAssembler<>::solution    solution;
//
//     geometryMap G = A.getMap(*m_mp);
//
//     space u = A.getSpace(dbasis);
//
//     A.initSystem();
//
//     gsFunctionExpr<> x("x",2);
//     gsFunctionExpr<> y("y",2);
//     variable ffx = A.getCoeff(x);
//     variable ffy = A.getCoeff(y);
//     auto j00 = fjac(ffx).tr()*jac(G)*fjac(ffx);
//     auto j10 = fjac(ffy).tr()*jac(G)*fjac(ffx);
//     auto j01 = fjac(ffx).tr()*jac(G)*fjac(ffy);
//     auto j11 = fjac(ffy).tr()*jac(G)*fjac(ffy);
//
//     auto g11 = j00*j00 + j10*j10;
//     auto g12 = j00*j01 + j10*j11;
//     auto g22 = j01*j01 + j11*j11;
//
//     auto uxi = grad(u)*fjac(ffx);
//     auto ueta = grad(u)*fjac(ffy);
//
//     auto detJ = j00*j11 - j01*j10;
//     auto detJinv = detJ.inv();
//     auto detJinv2 = detJinv*detJinv;
//         auto detJinv3 = detJinv2*detJinv;
//
//       auto d_detJ_dcx = uxi*j11 - ueta*j10 ;
//       auto d_g11pg22_dcx = 2*(uxi*j00 + ueta*j01);
//
//       auto fx = -d_detJ_dcx*(g11 + g22);
//       auto gx = d_g11pg22_dcx;
//
//       auto I = (g11+g22)*detJinv;
//
//       auto d_I_dcx = fx*detJinv2 + gx*detJinv;
//
//       auto d_detJinv2_dcx = -2 * d_detJ_dcx * detJinv3;
//       auto d_detJinv_dcx = - d_detJ_dcx * detJinv2;
//
//       auto d_fx_dcx_detJinv2 = -d_detJ_dcx*detJinv2*d_g11pg22_dcx.tr();
//       auto d_gx_dcx_detJinv = 2*(uxi*detJinv*uxi.tr() + ueta*detJinv*ueta.tr());
//
//       auto d2_I_dcx2 = fx * d_detJinv2_dcx.tr() + d_fx_dcx_detJinv2 + gx*d_detJinv_dcx.tr() + d_gx_dcx_detJinv;
//
//       A.assemble(-d2_I_dcx2);
//
//       gsMatrix<> xxMat = A.matrix();
//
//       A.initSystem();
//
//       auto d_detJ_dcy = ueta*j00 - uxi*j01 ;
//       auto d_g11pg22_dcy = 2*(uxi*j10 + ueta*j11);
//
//       auto fy = -d_detJ_dcy*(g11 + g22);
//       auto gy = d_g11pg22_dcy;
//
//
//       auto d_I_dcy = fy*detJinv2 + gy*detJinv;
//
//       auto d_detJinv2_dcy = -2 * d_detJ_dcy * detJinv3;
//       auto d_detJinv_dcy = - d_detJ_dcy * detJinv2;
//
//       auto d_fx_dcy_detJinv2 = ueta*(g11 + g22)*detJinv2*uxi.tr() - uxi*(g11 + g22)*detJinv2*ueta.tr() - d_detJ_dcx*detJinv2*d_g11pg22_dcy.tr();
//
//       auto d2_I_dcxcy = fx * d_detJinv2_dcy.tr() + d_fx_dcy_detJinv2 + gx*d_detJinv_dcy.tr();
//
//       A.assemble(-d2_I_dcxcy);
//       gsMatrix<> xyMat = A.matrix().transpose();
//
//       A.initSystem();
//
//       auto d_fy_dcy_detJinv2 = -d_detJ_dcy*detJinv2*d_g11pg22_dcy.tr();
//       auto d_gy_dcy_detJinv = 2*(uxi*detJinv*uxi.tr() + ueta*detJinv*ueta.tr());
//
//       auto d2_I_dcy2 = fy*d_detJinv2_dcy.tr() + d_fy_dcy_detJinv2 + gy*d_detJinv_dcy.tr() + d_gy_dcy_detJinv;
//
//       A.assemble(-d2_I_dcy2);
//
//
//       gsMatrix<> yyMat = A.matrix();
//
//         gsMatrix<> all(xxMat.rows() + yyMat.rows(),xxMat.cols() + yyMat.cols());
//
//         all << xxMat, xyMat.transpose(), xyMat, yyMat;
//
//         // Map twice, is there a better way?
//         gsMatrix<> all2 = mapMatrix(u.mapper(),all);
//
//         // Map tagged part
//         hessObjTagged = mapMatrixToTagged(u.mapper(),all2);
//
//         return mapMatrix(u.mapper(),all2);
//
//     }

// NEW IMPLENTATION FIXIT FINISH
    gsMatrix<> gsWinslow::hessObj(gsMatrix<> &hessObjTagged) const{
        gsExprAssembler<> A(1,1);
        gsMultiBasis<> dbasis(*m_mp);
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;

        geometryMap G = A.getMap(*m_mp);

        space u = A.getSpace(dbasis,m_mp->geoDim());

        A.initSystem();

        auto detJinv = jac(G).inv().det(); // The inverse of abs(det J)

        //A.assemble(2**detJinv
        //    - matrix_by_space(jac(G).inv(),jac(u)).trace()*
        //    (jac(G).tr()*jac(G)).trace().val()*detJinv);
        //
        A.assemble
        (
            2*(
                (jac(u)%jac(u).tr())*detJinv
                -   (jac(G).inv().tr()%jac(u))*(jac(G)%jac(u)).tr()*detJinv
                )
            +   (matrix_by_space(jac(G).inv(),jac(u))*(jac(G).inv()*jac(u))).trace()
                    *(jac(G)%jac(G))*detJinv
            -   2*matrix_by_space(jac(G).inv(),jac(u)).trace()*(jac(G)%jac(u))*detJinv
            +   matrix_by_space(jac(G).inv(),jac(u)).trace()
                    *matrix_by_space(jac(G).inv(),jac(u)).trace()
                    *(jac(G)%jac(G))*detJinv
        );

        gsMatrix<> all = A.matrix();

        // Map twice, is there a better way?
        gsMatrix<> all2 = mapMatrix(u.mapper(),all);

        // Map tagged part
        hessObjTagged = mapMatrixToTagged(u.mapper(),all2);

        return mapMatrix(u.mapper(),all2);

    }
