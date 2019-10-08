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

    gsMultiPatch<> x_mp = m_mp->coord(0);
    gsMultiPatch<> y_mp = m_mp->coord(1);

    variable x = A.getCoeff(x_mp);
    variable y = A.getCoeff(y_mp);

    auto Lx = (jac(G).tr()*jac(G)).adj() % hess(x);
    auto Ly = (jac(G).tr()*jac(G)).adj() % hess(y);

    real_t out = 0;

    out += ev.integral(lambda_1 * hess(G).sqNorm());
    out += ev.integral(lambda_2 * jac(G).sqNorm());
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

// Fixit, doesnt work yet.
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

    auto Lxd2Lxdcxcx = 2*(2*Lx.val()*dJacdc(u,0)%(dJacdc(u,0)*hess(x).adj())
                        + 2*Lx.val()*dJacdc(u,0)%(jac(G)*hess(u).adj()))
                        + 2*2*Lx.val()*hess(u)%(jac(G).tr()*dJacdc(u,0)).adj();

    auto Lxd2Lxdcxcy = 2*(2*Lx.val()*dJacdc(u,0)%(dJacdc(u,1)*hess(x).adj()))
                        + 2*(2*Lx.val()*jac(G).tr()*dJacdc(u,1)).adj()%hess(u);

    auto Lxd2Lxdcycy = 2*(2*Lx.val()*dJacdc(u,1)%(dJacdc(u,1)*hess(x).adj()));

    auto dLydcx = 2*(jac(G).tr()*dJacdc(u,0)).adj()%hess(y);
    auto dLydcy = 2*(jac(G).tr()*dJacdc(u,1)).adj()%hess(y) + hess(u) %(jac(G).tr()*jac(G)).adj();

    auto Lyd2Lydcycy = 2*(2*Ly.val()*dJacdc(u,1)%(dJacdc(u,1)*hess(y).adj())
                        + 2*Ly.val()*dJacdc(u,1)%(jac(G)*hess(u).adj()))
                        + 2*2*Ly.val()*hess(u)%(jac(G).tr()*dJacdc(u,1)).adj();

    auto Lyd2Lydcxcy = 2*(2*Ly.val()*dJacdc(u,0)%(dJacdc(u,1)*hess(y).adj()))
                        + 2*2*Ly.val()*hess(u)%(jac(G).tr()*dJacdc(u,0)).adj();

    auto Lyd2Lydcxcx = 2*(2*Ly.val()*dJacdc(u,0)%(dJacdc(u,0)*hess(y).adj()));

    // x x
    A.initSystem();

    A.assemble(2*lambda_1*hess(u)%hess(u));
    A.assemble(2*lambda_2*dJacdc(u,0)%dJacdc(u,0));
    A.assemble(2*dLxdcx*dLxdcx.tr());
    A.assemble(Lxd2Lxdcxcx);

    A.assemble(2*dLydcx*dLydcx.tr());
    A.assemble(Lyd2Lydcxcx);

    gsMatrix<> xxMat = A.matrix();
    gsInfo << "norm = " << xxMat.norm() << "\n";

    // x y
    A.initSystem();

    A.assemble(2*lambda_2*dJacdc(u,1)%dJacdc(u,0));
    A.assemble(2*dLxdcy*dLxdcx.tr());
    A.assemble(Lxd2Lxdcxcy);

    A.assemble(2*dLydcy*dLydcx.tr());
    A.assemble(Lyd2Lydcxcy);

    gsMatrix<> xyMat = A.matrix();
    gsInfo << "norm = " << xyMat.norm() << "\n";


    // y y
    A.initSystem();

    A.assemble(2*lambda_1*hess(u)%hess(u));
    A.assemble(2*lambda_2*dJacdc(u,1)%dJacdc(u,1));
    A.assemble(2*dLxdcy*dLxdcy.tr());
    A.assemble(Lxd2Lxdcycy);

    A.assemble(2*dLydcy*dLydcy.tr());
    A.assemble(Lyd2Lydcycy);

    gsMatrix<> yyMat = A.matrix();
    gsInfo << "norm = " << yyMat.norm() << "\n";


    gsMatrix<> all(xxMat.rows() + yyMat.rows(),xxMat.cols() + yyMat.cols());

    all << xxMat, xyMat.transpose(), xyMat, yyMat;

    // Map twice, is there a better way?
    gsMatrix<> all2 = mapMatrix(u.mapper(),all);

    // Map tagged part
    hessObjTagged = mapMatrixToTagged(u.mapper(),all2);

    return mapMatrix(u.mapper(),all2);

    }
