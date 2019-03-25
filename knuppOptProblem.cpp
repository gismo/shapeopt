#include <gismo.h>
#include <fstream>
#include "knuppOptProblem.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

real_t knuppOptProblem::evaluateOnPatch(index_t i) const{
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
  auto detJ2 = detJ*detJ;

  return ev.integral(k_A*detJ2 + k_O*g12*g12);
}

void knuppOptProblem::evaluateDerivOnPatch(index_t i, gsVector<> &xVec, gsVector<> &yVec) const{
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

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

  gsJacDetField<real_t> jacDetField(mp->patch(i));
  variable detJ = ev.getVariable(jacDetField);

  auto d_detJ_dcx = uxi*j11 - ueta*j10 ;
  auto d_detJ2_dcx = 2*d_detJ_dcx*detJ;

  auto d_g12_dcx = uxi*j01 + ueta*j00;
  auto d_g12m2_dcx = 2*d_g12_dcx*g12;

  A.assemble(k_A*d_detJ2_dcx + k_O*d_g12m2_dcx);

  xVec = A.rhs();

  A.initSystem();

  auto d_detJ_dcy = ueta*j00 - uxi*j01 ;
  auto d_detJ2_dcy = 2*d_detJ_dcy*detJ;

  auto d_g12_dcy = uxi*j11 + ueta*j10;
  auto d_g12m2_dcy = 2*d_g12_dcy*g12;

  A.assemble(k_A*d_detJ2_dcy + k_O*d_g12m2_dcy);
  yVec = A.rhs();


}

void knuppOptProblem::evaluate2ndDerivOnPatch(index_t i, gsMatrix<> &xxMat, gsMatrix<> &xyMat, gsMatrix<> &yyMat) const{
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

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

  gsJacDetField<real_t> jacDetField(mp->patch(i));
  variable detJ = ev.getVariable(jacDetField);

  auto d_detJ_dcx = uxi*j11 - ueta*j10 ;
  auto d2_detJ2_dcx2 = 2*d_detJ_dcx*d_detJ_dcx.tr();

  auto d_g12_dcx = uxi*j01 + ueta*j00;
  auto d2_g12_dcx2_g12 = uxi*g12*ueta.tr() + ueta*g12*uxi.tr();
  auto d2_g12m2_dcx2 = 2*(d_g12_dcx*d_g12_dcx.tr() + d2_g12_dcx2_g12);

  A.assemble(k_A*d2_detJ2_dcx2 + k_O*d2_g12m2_dcx2);

  xxMat = A.matrix();
  A.initSystem();

  auto d_detJ_dcy = ueta*j00 - uxi*j01 ;
  // auto d2_detJ_dcycx_detJ = ueta*detJ*uxi.tr() - uxi*detJ*ueta.tr() ;
  // auto d2_detJ2_dcycx = 2*(d_detJ_dcy*d_detJ_dcx.tr() + d2_detJ_dcycx_detJ) ;

  auto d_g12_dcy = uxi*j11 + ueta*j10;
  // auto d2_g12m2_dcycx = 2*d_g12_dcy*d_g12_dcx.tr();

  auto d2_detJ_dcxcy_detJ = uxi*detJ*ueta.tr() - ueta*detJ*uxi.tr();
  auto d2_detJ2_dcxcy = 2*(d_detJ_dcx*d_detJ_dcy.tr() + d2_detJ_dcxcy_detJ);

  auto d2_g12m2_dcxcy = 2*d_g12_dcx*d_g12_dcy.tr();
  A.assemble(k_A*d2_detJ2_dcxcy + k_O*d2_g12m2_dcxcy);

  xyMat = A.matrix().transpose();

  A.initSystem();

  auto d2_detJ2_dcy2 = 2*d_detJ_dcy*d_detJ_dcy.tr();

  auto d2_g12_dcy2_g12 = uxi*g12*ueta.tr() + ueta*g12*uxi.tr();
  auto d2_g12m2_dcy2 = 2*(d_g12_dcy*d_g12_dcy.tr() + d2_g12_dcy2_g12);

  A.assemble(k_A*d2_detJ2_dcy2 + k_O*d2_g12m2_dcy2);

  yyMat = A.matrix();

}
