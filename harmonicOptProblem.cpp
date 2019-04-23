#include <gismo.h>
#include <fstream>
#include "harmonicOptProblem.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

harmonicOptProblem::harmonicOptProblem(gsMultiPatch<>* mpin): paraWithMapper(mpin), dJC(mpin){
  setupMapper(); // Set mapper to fix degrees of freedoms and set current design
  m_numDesignVars = m_mappers[0].freeSize() + m_mappers[1].freeSize();
  m_curDesign = getDesignVariables();

  m_desLowerBounds.setOnes(m_numDesignVars);
  m_desLowerBounds *= -1e9;
  m_desUpperBounds.setOnes(m_numDesignVars);
  m_desUpperBounds *= 1e9;

  m_numConstraints = 0;

  computeJacStructure();
  m_conLowerBounds.setZero(m_numConstraints);
  m_conUpperBounds.setZero(m_numConstraints);
}

real_t harmonicOptProblem::evalObj() const{
  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(*mp);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  // ev.options().setReal("quA",quA);
  // ev.options().setInt("quB",quB);
  // gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(*mp);

  gsFunctionExpr<> f1("x", "-y",2);
  variable ffun = ev.getVariable(f1);
  auto D = fjac(ffun);

  // return ev.integral(jac(G).sqNorm());
  // return ev.integral(hess(G).sqNorm());
  return ev.integral((D*hess(G)*D*fform(G)).trace().sqNorm().val() + lambda_1 * hess(G).sqNorm() + lambda_2*jac(G).sqNorm());
}

gsVector<> harmonicOptProblem::gradientObj() const{
  // gsInfo << "evalObj\n";
  gsExprAssembler<> A(1,1);
  gsMultiBasis<> dbasis(*mp);
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  // ev.options().setReal("quA",quA);
  // ev.options().setInt("quB",quB);
  // gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  geometryMap G = A.getMap(*mp);

  space u = A.getSpace(dbasis); // Vector space for derivatives, for 2D

  gsFunctionExpr<> f1("x", "-y",2);
  variable ffun = ev.getVariable(f1);
  auto D = fjac(ffun);

  gsFunctionExpr<> fe1("1","0",2);
  variable e1 = ev.getVariable(fe1);

  A.initSystem();
  // A.assemble(2*matrix_by_space(jac(G).tr(),jac(u)).trace());
  A.assemble(2*e1.tr()*collapse(hess(G)[0].tr(),hess(u)));

  return A.rhs();
  // return ev.integral((D*hess(G)*D*fform(G)).trace().sqNorm().val() + lambda_1 * hess(G).sqNorm() + lambda_2*jac(G).sqNorm());
}

real_t harmonicOptProblem::evalObj( const gsAsConstVector<real_t> & u) const {
  // gsInfo << "...evalObj\n";
  updateDesignVariables(u);
  return evalObj();
}

void harmonicOptProblem::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
    result[0] = 0;
}

void harmonicOptProblem::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const {
}

void harmonicOptProblem::print(){
  gsInfo << "m_numDesignVars  = " <<  m_numDesignVars << "\n";
  gsInfo << "m_numConstraints  = " <<  m_numConstraints << "\n";

  gsInfo << "m_numConJacNonZero  = " <<  m_numConJacNonZero << "\n";

  gsInfo << "m_conUpperBounds.size() = " << m_conUpperBounds.size() << "\n";
  gsInfo << "m_conLowerBounds.size() = " << m_conLowerBounds.size() << "\n";
  gsInfo << "m_desUpperBounds.size() = " << m_desUpperBounds.size() << "\n";
  gsInfo << "m_desLowerBounds.size() = " << m_desLowerBounds.size() << "\n";

  gsInfo << "m_curDesign.size() = " << m_curDesign.size() << "\n";

  // gsMatrix<> disp(m_desUpperBounds.size(),2);
  // disp << m_desUpperBounds,m_desLowerBounds;
  // gsInfo << ".. design upper and lower bounds\n";
  // gsInfo << disp << "\n";
  //
  gsMatrix<> disp2(m_conUpperBounds.size(),2);
  disp2 << m_conUpperBounds,m_conLowerBounds;
  gsInfo << ".. constraint upper and lower bounds\n";
  gsInfo << disp2 << "\n";

}

// void paraOptProblem::writeToFile(gsVector<> vec, std::string name) const{
//     gsInfo << "WRITING to " << name << "\n";
//   std::ofstream f(name);
//   for (auto &e : vec) f << std::setprecision(12) << e << "\n";
// }
//
// void paraOptProblem::writeToFile(gsMatrix<> mat, std::string name) const{
//     gsInfo << "WRITING to " << name << "\n";
//   std::ofstream f(name);
//   for(index_t i = 0; i < mat.rows(); i++){
//     for(index_t j = 0; j < mat.cols(); j++){
//       f << std::setprecision(20) << mat(i,j) << " ";
//     }
//     f << "\n";
//   }
// }
//
// void paraOptProblem::loadFromFile(std::string name){
//   std::ifstream f(name);
//   gsVector<> des(m_numDesignVars);
//   for(index_t i = 0; i < m_numDesignVars; i++){
//       f >> des[i];
//   }
//   dJC.updateDesignVariables(des);
//
// }
//
// gsMatrix<> paraOptProblem::hessianObj() const{
//   gsMatrix<> hess;
//   hess.setZero(m_numDesignVars,m_numDesignVars);
//   index_t ind = 0;
//   for(index_t i = 0; i < mp->nBoxes(); i++){
//     // gsInfo << "\ni = " << i << "\n";
//     gsMatrix<> xxMat,xyMat,yyMat;
//     evaluate2ndDerivOnPatch(i,xxMat,xyMat,yyMat);
//     index_t ncoefs = mp->patch(i).coefsSize();
//     hess.block(ind,ind,ncoefs,ncoefs) = xxMat;
//     hess.block(ind,ind+dJC.n_controlpoints,ncoefs,ncoefs) = xyMat.transpose();
//     hess.block(ind+dJC.n_controlpoints,ind,ncoefs,ncoefs) = xyMat;
//     hess.block(ind+dJC.n_controlpoints,ind+dJC.n_controlpoints,ncoefs,ncoefs) = yyMat;
//     ind += ncoefs;
//   }
//   return hess;
// }
//
// gsVector<> paraOptProblem::gradientObj() const{
//   gsVector<> grad;
//   grad.setZero(m_numDesignVars);
//   index_t ind = 0;
//   for(index_t i = 0; i < mp->nBoxes(); i++){
//     // gsInfo << "\ni = " << i << "\n";
//     gsVector<> xVec,yVec;
//     evaluateDerivOnPatch(i,xVec,yVec);
//     index_t ncoefs = mp->patch(i).coefsSize();
//     grad.segment(ind,ncoefs) = xVec;
//     grad.segment(ind+dJC.n_controlpoints,ncoefs) = yVec;
//     ind += ncoefs;
//   }
//   return grad;
// }

void harmonicOptProblem::reset(){
  // Mapper is already setup
  m_curDesign = getDesignVariables();
}
