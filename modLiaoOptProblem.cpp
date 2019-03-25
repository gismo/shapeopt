#include <gismo.h>
#include <fstream>
#include "modLiaoOptProblem.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;


real_t modLiaoOptProblem::evaluateOnPatch(index_t i) const{
  real_t result = 0;
  gsMultiPatch<> singlePatch(mp->patch(i));

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(singlePatch);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;
  
  geometryMap G = A.getMap(singlePatch);

  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = ev.getVariable(x);
  variable fy = ev.getVariable(y);
  auto j00 = grad(fx)*jac(G)*grad(fx).tr();
  auto j10 = grad(fy)*jac(G)*grad(fx).tr();
  auto j01 = grad(fx)*jac(G)*grad(fy).tr();
  auto j11 = grad(fy)*jac(G)*grad(fy).tr();

  auto g11 = j00*j00 + j10*j10;
  auto g12 = j00*j01 + j10*j11;
  auto g22 = j01*j01 + j11*j11;

  gsJacDetField<real_t> jacDetField(mp->patch(i));
  variable detJ = ev.getVariable(jacDetField);
  auto detJinv = detJ.inv();

  return ev.integral(detJinv*detJinv*(g11 + g22)*(g11 + g22));
}

void modLiaoOptProblem::evaluateDerivOnPatch(index_t i, gsVector<> &xVec, gsVector<> &yVec) const{
  gsMultiPatch<> singlePatch(mp->patch(i));

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(singlePatch);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(singlePatch);

  space u = A.getSpace(dbasis);

  A.initSystem();

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

  gsJacDetField<real_t> jacDetField(mp->patch(i));
  variable detJ = ev.getVariable(jacDetField);
  auto detJinv = detJ.inv();
  auto detJinv2 = detJinv*detJinv;

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

  auto dIdcx = (ueta*j10 - uxi*j11)*(g11 + g22)*detJinv2 + 2*(uxi*j00 + ueta*j01)*detJinv;
  A.assemble(2*dIdcx*(g11+g22)*detJinv);

  xVec = A.rhs();

  A.initSystem();

  auto dIdcy = (uxi*j01 - ueta*j00)*(g11 + g22)*detJinv2 + 2*(uxi*j10 + ueta*j11)*detJinv;
  A.assemble(2*dIdcy*(g11+g22)*detJinv);
  yVec = A.rhs();


}

void modLiaoOptProblem::evaluate2ndDerivOnPatch(index_t i, gsMatrix<> &xxMat, gsMatrix<> &xyMat, gsMatrix<> &yyMat) const{
  gsMultiPatch<> singlePatch(mp->patch(i));

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(singlePatch);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(singlePatch);

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

  gsJacDetField<real_t> jacDetField(mp->patch(i));
  variable detJ = ev.getVariable(jacDetField);
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

  auto d_fx_dcx_detJinv2_I = -d_detJ_dcx*detJinv2*I*d_g11pg22_dcx.tr();
  auto d_gx_dcx_detJinv_I = 2*(uxi*detJinv*I*uxi.tr() + ueta*detJinv*I*ueta.tr());

  auto d2_I_dcx2_I = fx * I * d_detJinv2_dcx.tr() + d_fx_dcx_detJinv2_I + gx*I*d_detJinv_dcx.tr() + d_gx_dcx_detJinv_I;

  A.assemble(2*(d2_I_dcx2_I + d_I_dcx*d_I_dcx.tr()));

  xxMat = A.matrix();

  A.initSystem();

  auto d_detJ_dcy = ueta*j00 - uxi*j01 ;
  auto d_g11pg22_dcy = 2*(uxi*j10 + ueta*j11);

  auto fy = -d_detJ_dcy*(g11 + g22);
  auto gy = d_g11pg22_dcy;


  auto d_I_dcy = fy*detJinv2 + gy*detJinv;

  auto d_detJinv2_dcy = -2 * d_detJ_dcy * detJinv3;
  auto d_detJinv_dcy = - d_detJ_dcy * detJinv2;

  auto d_fx_dcy_detJinv2_I = ueta*(g11 + g22)*detJinv2*I*uxi.tr() - uxi*(g11 + g22)*detJinv2*I*ueta.tr() - d_detJ_dcx*detJinv2*I*d_g11pg22_dcy.tr();

  auto d2_I_dcxcy_I = fx * I * d_detJinv2_dcy.tr() + d_fx_dcy_detJinv2_I + gx*I*d_detJinv_dcy.tr();

  A.assemble(2*(d2_I_dcxcy_I + d_I_dcx*d_I_dcy.tr()));
  xyMat = A.matrix().transpose();

  A.initSystem();

  auto d_fy_dcy_detJinv2_I = -d_detJ_dcy*detJinv2*I*d_g11pg22_dcy.tr();
  auto d_gy_dcy_detJinv_I = 2*(uxi*detJinv*I*uxi.tr() + ueta*detJinv*I*ueta.tr());

  auto d2_I_dcy2_I = fy*I*d_detJinv2_dcy.tr() + d_fy_dcy_detJinv2_I + gy*I*d_detJinv_dcy.tr() + d_gy_dcy_detJinv_I;

  A.assemble(2*(d2_I_dcy2_I + d_I_dcy*d_I_dcy.tr()));


  yyMat = A.matrix();

}
