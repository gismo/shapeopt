#include <gismo.h>
#include <fstream>
#include "liaoOptProblem.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

real_t liaoOptProblem::evaluateOnPatch(index_t i) const{
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
  return ev.integral(g11*g11 + 2*g12*g12 + g22*g22);
}

void liaoOptProblem::evaluateDerivOnPatch(index_t i, gsVector<> &xVec, gsVector<> &yVec) const{
  gsMultiPatch<> singlePatch(mp->patch(i));

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(singlePatch);
  A.setIntegrationElements(dbasis);

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

  auto h1x = j00*g11 + j01*g12;
  auto h2x = j01*g22 + j00*g12;

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

  A.assemble(4*uxi*h1x + 4*ueta*h2x);

  xVec = A.rhs();

  A.initSystem();

  auto h1y = j10*g11 + j11*g12;
  auto h2y = j11*g22 + j10*g12;

  A.assemble(4*uxi *h1y+ 4*ueta*h2y);
  yVec = A.rhs();


}

void liaoOptProblem::evaluate2ndDerivOnPatch(index_t i, gsMatrix<> &xxMat, gsMatrix<> &xyMat, gsMatrix<> &yyMat) const{
  gsMultiPatch<> singlePatch(mp->patch(i));

  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(singlePatch);
  A.setIntegrationElements(dbasis);

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

  auto uxi = grad(u)*fjac(fx);
  auto ueta = grad(u)*fjac(fy);

  auto dh1xdx = uxi*(3*j00*j00 + j10*j10 + j01*j01) + ueta*(2*j00*j01 + j10*j11);
  auto dh2xdx = ueta*(3*j01*j01 + j11*j11 + j00*j00) + uxi*(2*j00*j01 + j10*j11);

  A.assemble(4*dh1xdx*uxi.tr() + 4*dh2xdx*ueta.tr());

  xxMat = A.matrix();

  A.initSystem();

  auto dh1xdy = uxi*(2*j00*j10 + j01*j11) + ueta*j01*j10;
  auto dh2xdy = ueta*(j00*j10 + 2*j01*j11) + uxi*j11*j00;

  A.assemble(4*dh1xdy*uxi.tr() + 4*dh2xdy*ueta.tr());
  xyMat = A.matrix();

  A.initSystem();

  auto dh1ydy = uxi*(3*j10*j10 + j00*j00 + j11*j11) + ueta*(2*j10*j11 + j00*j01);
  auto dh2ydy = ueta*(3*j11*j11 + j01*j01 + j10*j10) + uxi*(2*j10*j11 + j00*j01);

  A.assemble(4*dh1ydy*uxi.tr() + 4*dh2ydy*ueta.tr());

  yyMat = A.matrix();

}
