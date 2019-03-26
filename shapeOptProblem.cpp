#include <gismo.h>
#include "shapeOptProblem.h"

shapeOptProblem::shapeOptProblem(gsMultiPatch<>* mpin): mp(mpin), dJC(mpin), iC(mpin), SE(mpin), pOP(mpin), linOP(&pOP){
  // Calculate number of design variables OBS assuming same n.o. controlpoints in each direction
  index_t n_cc = mp->patch(antennaPatch).coefsSize();
  index_t n_coefsPs = sqrt(n_cc);
  m_numDesignVars = n_coefsPs*4*2; //4 sides on patch, 2 means both x and y

  // Calculate number of constraints
  m_numConstraints = dJC.n_constraints;

  m_conLowerBounds.setOnes(m_numConstraints);
  m_conUpperBounds.setOnes(m_numConstraints);

  m_conUpperBounds.segment(0,dJC.n_constraints) = dJC.getUpperBounds(m_eps);
  m_conLowerBounds.segment(0,dJC.n_constraints) *= -iC.aBigNumber;

  // Find design indicies, go through all interfaces and check if bnd to patch
  // FIXIT: you can get patchSides of geometry by copying it to a gsMultiPatch obj.
  designIndiciesLocal.setZero(m_numDesignVars/2);
  designIndiciesGlobal.setZero(m_numDesignVars);
  index_t j = 0;
  for(index_t i = 0; i < mp->nInterfaces(); i++){
    boundaryInterface bI = mp->bInterface(i);
    patchSide ps;
    if (bI.first().patch == antennaPatch){
      ps = bI.first();
    } else if (bI.second().patch == antennaPatch){
      ps = bI.second();
    } else {
      continue;
    }
    gsVector<index_t> psIndicies(n_coefsPs);
    iC.getLocalPatchSideIndicies(ps,psIndicies);
    designIndiciesLocal.segment(j,n_coefsPs) = psIndicies;
    designIndiciesGlobal.segment(j,n_coefsPs) = iC.getPatchSideIndicies(ps);
    j += n_coefsPs;
  }

  for(index_t i = 0; i < m_numDesignVars/2;i++){
    designIndiciesGlobal[i+m_numDesignVars/2] = designIndiciesGlobal[i] + dJC.n_controlpoints;
  }

  m_curDesign.setZero(m_numDesignVars,1);

  gsVector<> xRef(m_numDesignVars);
  gsVector<> des = pOP.getDesignVariables();
  for(index_t i = 0; i < m_numDesignVars; i++){
    xRef[i] = des[designIndiciesGlobal[i]];
  }

  // gsInfo << "xRef:\n" << xRef << "\n";

  m_desLowerBounds.setOnes(m_numDesignVars);
  m_desUpperBounds.setOnes(m_numDesignVars);

  real_t lb_x = -d_bw*SE.pde_L_f/2;
  real_t ub_x = d_bw*SE.pde_L_f/2;
  real_t lb_y = d_g/2*SE.pde_L_f;
  real_t ub_y = (d_g/2 + d_bh)*SE.pde_L_f;

  gsInfo << "\n\nlbx: " << lb_x;
  gsInfo << "\nubx: " << ub_x;
  gsInfo << "\nlby: " << lb_y;
  gsInfo << "\nuby: " << ub_y << "\n";

  for(index_t i = 0; i < m_numDesignVars/2; i++){
    m_desLowerBounds[i] = lb_x - xRef[i];
    m_desUpperBounds[i] = ub_x - xRef[i];

    m_desLowerBounds[i+m_numDesignVars/2] = lb_y - xRef[i+m_numDesignVars/2];
    m_desUpperBounds[i+m_numDesignVars/2] = ub_y - xRef[i+m_numDesignVars/2];
  }

  // gsInfo << "desLowerBounds:\n" << m_desLowerBounds << "\n";
  // gsInfo << "desUpperBounds:\n" << m_desUpperBounds << "\n";
  // gsVector<> des;
  // des.setZero(m_numDesignVars);
  // for (index_t i = 0; i < m_numDesignVars; i++){
  //   des[i] += 1;
  // }
  // gsVector<> gradCon = derivVolumeOfPatch(antennaPatch);
  // gsMatrix<> gradTrans = gradCon.transpose();
  // IpOptSparseMatrix J(gradTrans,-1); // Jacobiant of volume constraint

  // Call getDvectors to generate LU factorization and stuff
  gsVector<> tmp(dJC.n_constraints);
  dJC.getDvectors(tmp);
  // gsInfo << "D vectors :" << tmp << "\n";
  IpOptSparseMatrix J = derivativeOfDetJac(); // Jacobiant of detJ constraints
  // J.concatenate(J2,"col");

  m_numConJacNonZero = J.nnz();
  m_conJacRows = J.rows();
  m_conJacCols = J.cols();

  gsInfo << "m_numConstraints " << m_numConstraints << "\n";
  gsInfo << "m_numConJacNonZero " << m_numConJacNonZero << "\n";

  // Setup obj fun
  delta = *(new gsFunctionExpr<>("exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));
  ddeltadx = *(new gsFunctionExpr<>("-x/(0.1^2)*exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));
  ddeltady = *(new gsFunctionExpr<>("-y/(0.1^2)*exp(-x^2/(2*0.1^2) -y^2/(2*0.1^2))",2));

  gsInfo << "\n\nevalObj() : " << evalObj() << "\n";
  gsInfo << "volumeOfPatch() : " << volumeOfPatch(antennaPatch) << "\n\n";
  // writeToFile(dJC.getDesignVariables(),"all_cps_before.txt");
  // updateDesignVariables(des);
  // writeToFile(dJC.getDesignVariables(),"all_cps_after.txt");

}

real_t shapeOptProblem::evalObj() const {
  gsMultiPatch<> u_real,u_imag;
  SE.solve(u_real,u_imag);

  gsExprAssembler<> A(1,1);
  gsExprEvaluator<> ev(A);

  gsOptionList opts = A.options();
  opts.setInt("quB",20);
  opts.setReal("quA",20);
  A.setOptions(opts);

  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(*mp);
  gsMultiBasis<> dbasis(SE.dbasis);
  A.setIntegrationElements(dbasis);

  // variable deltaf = A.getCoeff(delta,G);
  variable df = A.getCoeff(delta,G);

  // gsInfo<<"Plotting in Paraview...\n";
  // ev.options().setSwitch("plot.elements", true);
  // ev.writeParaview( deltaf    , G, "deltaf");
  variable u_r = A.getCoeff(u_real);
  variable u_i = A.getCoeff(u_imag);

  // return -2*ev.integral(df*meas(G));
  // return -2*ev.integral(df*u_r*meas(G));
  return -2*ev.integral(df*(u_r.sqNorm() + u_i.sqNorm())*meas(G));

}

real_t shapeOptProblem::evalObj( const gsAsConstVector<real_t> & u ) const {
  // gsInfo << "evalObj\n";
  updateDesignVariables(u);
  return evalObj();
  // return 0;
}

gsVector<> shapeOptProblem::gradientObjWithoutAdjoint() const{
  gsMultiPatch<> u_real,u_imag;
  SE.solve(u_real,u_imag);

  gsMatrix<> dcdx = derivativeOfDesignUpdate();
  gsMatrix<> dEdc = -2*(evaluateDerivativeTerm1(u_real,u_imag) + evaluateDerivativeTerm2(u_real,u_imag));

  gsMatrix<> out = dEdc.transpose()*dcdx;

  // gsInfo << "(" << out.rows() <<", " << out.cols() << ")\n" << std::flush;
  return out.transpose();
};

// Adjoint method for sensitivities.
gsVector<> shapeOptProblem::gradientObj() const{
  gsMultiPatch<> u_real,u_imag;
  SE.solve(u_real,u_imag);

  gsMatrix<> dcdx = derivativeOfDesignUpdate();

  gsMatrix<> mat = SE.getDerivativeWithoutSolving();
  gsVector<> dJdu = getObjDerivativeDu(u_real,u_imag);

  gsVector<> adjoint = SE.solveAdjoint(dJdu);

  // SE.printMatSize(mat,"mat");
  // SE.printMatSize(adjoint,"adjoint");
  gsVector<> term2 = adjoint.transpose()*mat;

  gsMatrix<> dEdc = -2*(evaluateDerivativeTerm1(u_real,u_imag) + term2);

  gsMatrix<> out = dEdc.transpose()*dcdx;

  // gsInfo << "(" << out.rows() <<", " << out.cols() << ")\n" << std::flush;
  return out.transpose();
};

gsVector<> shapeOptProblem::evaluateDerivativeTerm1(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(*mp);
  A.setIntegrationElements(SE.dbasis);
  gsExprEvaluator<> ev(A);

  // opts.setReal("quA",2);
  gsOptionList opts = A.options();
  opts.setInt("quB",1);
  // gsInfo << opts << "\n";
  A.setOptions(opts);

  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(*mp);

  space u = A.getSpace(dbasis);

  variable u_r = A.getCoeff(u_real);
  variable u_i = A.getCoeff(u_imag);

  variable df = A.getCoeff(delta,G);
  variable ddfdx = A.getCoeff(ddeltadx,G);
  variable ddfdy = A.getCoeff(ddeltady,G);

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

gsVector<> shapeOptProblem::evaluateDerivativeTerm2(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
  // Perhaps the solve is called unecessary in here: FIXIT
  gsMatrix<> dudc = SE.getDerivativeOfU();

  gsVector<> dJdu = getObjDerivativeDu(u_real,u_imag);

  // gsInfo << "\n\n-------RETURN SOME STUFF-------\n\n";
  // SE.printMatSize(dudc,"dudc");
  // SE.printMatSize(dJdu,"dJdu");
  return dudc*dJdu;
}

gsVector<> shapeOptProblem::getObjDerivativeDu(gsMultiPatch<> &u_real, gsMultiPatch<> &u_imag) const{
  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(*mp);
  A.setIntegrationElements(SE.dbasis);
  gsExprEvaluator<> ev(A);


  // opts.setReal("quA",2);
  gsOptionList opts = A.options();
  opts.setInt("quB",1);
  // gsInfo << opts << "\n";
  A.setOptions(opts);

  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(*mp);

  space du = A.getSpace(SE.dbasis);
  du.setInterfaceCont(0);

  variable u_r = A.getCoeff(u_real);
  variable u_i = A.getCoeff(u_imag);

  variable df = A.getCoeff(delta,G);

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

gsVector<> shapeOptProblem::evalCon() const {
  gsVector<> tmp(dJC.n_constraints);
  dJC.getDvectors(tmp);

  // gsInfo << "\n...Max of D vector : " << tmp.maxCoeff() << "\n";

  return tmp;

}

void shapeOptProblem::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  updateDesignVariables(u);
  // gsInfo << "evalCon_into\n" << std::flush;
  // gsInfo << "...evalCon_into" << "\n";

  result = evalCon();
}

void shapeOptProblem::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
  updateDesignVariables(u);
  // gsInfo << "jacobCon_into\n";
  // gsInfo << "results size : (" << result.rows() << ", " << result.cols() << ")\n";
  // gsVector<> gradCon = derivVolumeOfPatch(antennaPatch);
  // gsMatrix<> gradTrans = gradCon.transpose();
  // IpOptSparseMatrix J(gradTrans,-1); // Jacobiant of volume constraint
  char str [50];

  if (counter1 >= 0){
    sprintf(str,"../results/shapeopt1/design_%d.txt",counter1++);
    writeToFile(dJC.getDesignVariables(),std::string(str));

    // sprintf(str,"shapeOptProblemGradTest12/x_%d.txt",counter1);
    // writeToFile(u,std::string(str));
  } else {
    counter1++;
  }

  IpOptSparseMatrix J = derivativeOfDetJac(); // Jacobiant of detJ constraints

  // J.concatenate(J2,"col");
  // gsInfo << "\nJ.nnz() : " << J.nnz() << "\n";
  // gsInfo << "J.values() size : (" << J.values().rows() << ", " << J.values().cols() << ")\n";
  result = J.values();
}

gsVector<> shapeOptProblem::getDesignVariables() const{
  return m_curDesign;
}

void shapeOptProblem::setCurrentDesign(gsVector<> x){
  m_curDesign = x;
};

void shapeOptProblem::updateDesignVariables(gsVector<> u) const {
  gsVector<> deltaCps;
  deltaCps.setZero(linOP.numDesignVars());
  for (index_t j = 0; j < m_numDesignVars; j++){
    deltaCps(designIndiciesGlobal[j]) = u[j];
  }
  linOP.solveAndUpdate(deltaCps);
}

gsVector<> shapeOptProblem::getUpdateToCps(gsVector<> u) const {
  gsVector<> deltaCps;
  deltaCps.setZero(linOP.numDesignVars());
  for (index_t j = 0; j < m_numDesignVars; j++){
    deltaCps(designIndiciesGlobal[j]) = u[j];
  }
  return linOP.solve(deltaCps);
}

void shapeOptProblem::writeToFile(gsVector<> vec, std::string name) const{
  std::ofstream f(name);
  for (auto &e : vec) f << std::setprecision(12) << e << "\n";
}

real_t shapeOptProblem::volumeOfPatch(index_t p) const {
  gsExprAssembler<> A(1,1);

  gsExprEvaluator<> ev(A);

  gsMultiPatch<> singlePatch(mp->patch(p));

  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(singlePatch);
  gsMultiBasis<> dbasis(singlePatch);
  A.setIntegrationElements(dbasis);

  return ev.integral(1.0*meas(G));
}

gsMatrix<> shapeOptProblem::derivVolumeOfPatch(index_t p) const {
  gsExprAssembler<> A(1,1);

  gsExprEvaluator<> ev(A);

  gsMultiPatch<> singlePatch(mp->patch(p));

  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(singlePatch);
  gsMultiBasis<> dbasis(singlePatch);
  A.setIntegrationElements(dbasis);

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

  auto detJ = j00*j11 - j10*j01;
	auto detJinv = 1/detJ.val();
	auto signOfDetJ = detJinv*meas(G);

  auto d_detJ_dcx = signOfDetJ*(uxi*j11 - ueta*j10) ;
  auto d_detJ_dcy = signOfDetJ*(ueta*j00 - uxi*j01) ;

  A.assemble(d_detJ_dcx);

  gsVector<> xVec = A.rhs();

  A.initSystem();

  A.assemble(d_detJ_dcy);

  gsVector<> yVec = A.rhs();

  gsMatrix<> dfdc(xVec.rows()*2,1);
  dfdc << xVec, yVec;

  gsMatrix<> tmp = derivativeOfDesignUpdate();
  gsMatrix<> dcdx = extractPatchRows(antennaPatch,tmp);

  // gsInfo << "(" << dfdc.rows() << ", " << dfdc.cols() << ")\n";
  // gsInfo << "(" << dcdx.rows() << ", " << dcdx.cols() << ")\n" << std::flush;

  gsMatrix<> out = dfdc.transpose()*dcdx;
  // gsInfo << "(" << out.rows() <<", " << out.cols() << ")\n" << std::flush;
  return out.transpose();
}

gsMatrix<> shapeOptProblem::derivativeOfDesignUpdate() const {
  gsMatrix<> out(linOP.numDesignVars(),m_numDesignVars); //the size is n.o. designvars for linearizedOptProblem times n.o. design variables for this problem

  gsVector<> zers;
  zers.setZero(numDesignVars());
  gsVector<> c_zero = getUpdateToCps(zers);
  for(index_t i = 0; i < m_numDesignVars; i++){
    gsVector<> ei;
    ei.setZero(m_numDesignVars);
    ei[i] = 1;
    gsVector<> ci = getUpdateToCps(ei) - c_zero; // Substract value at zero since it is an affine function
    for(index_t j = 0; j < linOP.numDesignVars(); j++){
      out(j,i) = ci[j];
    }
  }

  return out;
}

gsMatrix<> shapeOptProblem::extractPatchRows(index_t p, gsMatrix<> in) const {
  index_t ncc = mp->patch(p).coefsSize();
  index_t patchOffset = 0;
  for(int i = 1; i <= p; i++){
    patchOffset += mp->patch(i).coefsSize();
  }

  gsMatrix<> out(ncc*2,in.cols());

  out.block(0,0,ncc,in.cols()) = in.block(patchOffset,0,ncc,in.cols());
  out.block(ncc,0,ncc,in.cols()) = in.block(patchOffset + dJC.n_controlpoints,0,ncc,in.cols());
  return out;
}

IpOptSparseMatrix shapeOptProblem::derivativeOfDetJac() const{
  // FIXIT sparse structure of dDdc is lost.. Recover possibly?
  gsMatrix<> dDdc = dJC.getJacobian().asDense();
  gsMatrix<> dcdx = derivativeOfDesignUpdate();

  gsMatrix<> dDdx = dDdc*dcdx;
  IpOptSparseMatrix out(dDdx,-1);
  return out;
}

void shapeOptProblem::resetParametrizationToReference(){
  pOP.updateDesignVariables(linOP.refControlPoints());
}

void shapeOptProblem::gradObj_into ( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result ) const{
  updateDesignVariables(u);
  result = gradientObj();
}
