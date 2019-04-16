#include <gismo.h>
#include "detJacConstraint.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

detJacConstraint::detJacConstraint(gsMultiPatch<>* mpin): mp(mpin), m_detJacBasis(*mp),
        solversMassMatrix(mpin->nBoxes()),m_areSolversSetup(mpin->nBoxes()){
    for (index_t i = 0; i < mp->nBoxes(); i++){
        m_areSolversSetup[i] = false;
    }
    // Prepare basis for detJac
    // by setting the degree to (2p-1)
    int p = m_detJacBasis.maxCwiseDegree();
    m_detJacBasis.setDegree(2*p-1);
    // by reducing the continuity
    m_detJacBasis.reduceContinuity(1);
    // gsInfo << "\n..detJacBasis.size: " << m_detJacBasis.size(0) << "\n";
    // gsInfo << "\n..detJacBasis degree is: " << 2*p-1 << "\n";

    // Count the total number of controlpoints
    n_controlpoints = 0;
    n_constraints = 0;
    for(int i = 0; i < mp->nBoxes(); i++){
        n_controlpoints += mp->patch(i).coefsSize();
        n_constraints += m_detJacBasis.size(i);
        // gsInfo << "Patch " << i << " has " << mp->patch(i).coefsSize() << " cps... \n";
        // gsInfo << "Patch " << i << " gets " << m_detJacBasis.size(i) << " constraints... \n";
    }

    // gsInfo << "n_controlpoints = " << n_controlpoints << "\n \n";
    // gsInfo << "n_constraints = " << n_controlpoints << "\n \n";
}

gsVector<> detJacConstraint::generateDResultVector(){
    int len = 0;
    for(int i = 0; i < mp->nBoxes(); i++){
        len += m_detJacBasis.size(i);
    }
    gsVector<> result(len);
    return result;
}

void detJacConstraint::getDvectors(gsVector<> &result){
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    index_t start = 0;
    for(int i = 0; i < mp->nBoxes(); i++){
        gsExprAssembler<> A(1,1);

        gsMultiBasis<> dbasis(m_detJacBasis.basis(i));
        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);

        space u = A.getSpace(dbasis);

        A.initSystem();
        if (! m_areSolversSetup[i]){
            A.assemble(u*u.tr());
            solversMassMatrix[i].compute(A.matrix());
            m_areSolversSetup[i] = true;
        }

        gsMatrix<> solVector;
        solution u_sol = A.getSolution(u, solVector);

        // gsJacDetField<real_t> jacDetField(mp->patch(i));
        // variable detJ = ev.getVariable(jacDetField);
        // A.assemble(u*detJ);

        geometryMap G = A.getMap(mp->patch(i));
        A.assemble(u*jac(G).det());

        solVector = solversMassMatrix[i].solve(A.rhs());

        // Save result in result vector
        index_t len   = m_detJacBasis.size(i);
        result.segment(start,len) = solVector;

        start += m_detJacBasis.size(i);
    }



}

gsVector<> detJacConstraint::getDvectors(){
    gsVector<> out = generateDResultVector();
    getDvectors(out);
    return out;

}

void detJacConstraint::getDerivRhsFromPatch(index_t patch, gsSparseMatrix<> &xJac, gsSparseMatrix<> &yJac){
    gsMultiPatch<> singlePatch(mp->patch(patch));

    gsMultiBasis<> dJbas(m_detJacBasis.basis(patch));

    // Prepare assembler
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(singlePatch);
    A.setIntegrationElements(dJbas);
    gsExprEvaluator<> ev(A);

    // Define types
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(singlePatch);

    // Setup and assemble the two matricies
    space u = A.getSpace(dbasis,mp->geoDim());
    space v = A.getTestSpace(u,dJbas,1);
    // space v = A.getSpace(dJbas,1);

    // gsInfo << "dJbas.size : " << dJbas.size() << "\n";
    // gsInfo << "dbasis.size : " << dbasis.size() << "\n";

    A.initSystem();
    A.assemble(v*matrix_by_space(jac(G).inv(),jac(u)).trace().tr()*jac(G).det());
    // A.assemble(matrix_by_space(jac(G).inv(),jac(u)).trace()*jac(G).det());
    // A.assemble(v);
    // A.assemble(v*u.tr());

    // gsInfo << "\nsize (" << A.rhs().rows() << ", " << A.rhs().cols() << ")\n";
    // gsInfo << "\nsize (" << A.matrix().rows() << ", " << A.matrix().cols() << ")\n";
    index_t r = A.matrix().rows();
    index_t c = A.matrix().cols();
    xJac = A.matrix().block(0,0,r,c/2);
    yJac = A.matrix().block(0,c/2,r,c/2);

}

void detJacConstraint::getJacobianFromPatch(index_t patch, gsMatrix<> &xJac, gsMatrix<> &yJac){
    GISMO_ASSERT(m_areSolversSetup[patch],"Solver is not setup before calling detJacConstraint::getJacobianFromPatch");

    gsSparseMatrix<> xDrhs,yDrhs;
    getDerivRhsFromPatch(patch,xDrhs,yDrhs);

    xJac = solversMassMatrix[patch].solve(xDrhs);

    yJac = solversMassMatrix[patch].solve(yDrhs);
}

IpOptSparseMatrix detJacConstraint::getJacobian(){
    // For each patch generate Jacobian with respect to x and y coordinates of geometry
    memory::unique_ptr<IpOptSparseMatrix> xMat,yMat;         // Store a pointer to the object, to avoid calling the constructor
    for(index_t i = 0; i < mp->nBoxes(); i++){
        // get xJac and yJac for patch i
        gsMatrix<> xJac,yJac;
        getJacobianFromPatch(i,xJac,yJac);
        if (i == 0){
            xMat.reset(new IpOptSparseMatrix(xJac,-1));    // -1 indicates that IpOptSparseMatrix should generate dense matrix
            yMat.reset(new IpOptSparseMatrix(yJac,-1));    // --||-- ...
        } else {
            xMat->concatenate(IpOptSparseMatrix(xJac,-1),"diag");
            yMat->concatenate(IpOptSparseMatrix(yJac,-1),"diag");
        }
    }
    xMat->concatenate(*yMat,"row");       // Store full jacobian in xMat

    return *xMat;

}

gsVector<> detJacConstraint::getDesignVariables(){
    gsVector<> cx(n_controlpoints);
    gsVector<> cy(n_controlpoints);
    index_t j = 0;
    for(index_t i = 0; i < mp->nBoxes(); i++){
        for(index_t k = 0; k < mp->patch(i).coefsSize(); k++){
            cx[j] = mp->patch(i).coef(k,0);
            cy[j] = mp->patch(i).coef(k,1);
            j++;
        }
    }

    gsVector<> out(2*n_controlpoints);
    out << cx,cy;

    return out;
}

void detJacConstraint::updateDesignVariables(gsVector<> des){

    GISMO_ASSERT(des.size() == 2*n_controlpoints, "Design vector is of wrong size.!");

    gsVector<> cx = des.segment(0,n_controlpoints);
    gsVector<> cy = des.segment(n_controlpoints,n_controlpoints);

    index_t seg = 0;
    for(index_t i = 0; i < mp->nBoxes(); i++){
        index_t n_coefs = mp->patch(i).coefsSize();
        gsMatrix<> cc(n_coefs,2);

        cc << cx.segment(seg,n_coefs),
        cy.segment(seg,n_coefs);

        seg += n_coefs;

        mp->patch(i).setCoefs(cc);
    }
}

//FIXIT: detJacConstraint is coded for negative determinant.. FIXIT
gsVector<> detJacConstraint::getUpperBounds(real_t eps){
    gsVector<> out;
    out.setOnes(n_constraints);
    out *= -eps;
    return out;
}

void detJacConstraint::plotDetJ(std::string name){
  // gsInfo << "getDvectors\n" << std::flush;
  //gsInfo<<"Active options:\n"<< A.options() <<"\n";
  typedef gsExprAssembler<>::geometryMap geometryMap;
  typedef gsExprAssembler<>::variable    variable;
  typedef gsExprAssembler<>::space       space;
  typedef gsExprAssembler<>::solution    solution;

  gsExprAssembler<> A(1,1);

  // Elements used for numerical integration
  A.setIntegrationElements(m_detJacBasis);
  gsExprEvaluator<> ev(A);

  space u = A.getSpace(m_detJacBasis);

  geometryMap G = A.getMap(*mp);

  A.initSystem();
  A.assemble(u*u.tr());
  gsSparseSolver<>::CGDiagonal solver;
  solver.compute(A.matrix());

  gsMatrix<> solVector;
  solution u_sol = A.getSolution(u, solVector);

  gsFunctionExpr<> x("x",2);
  gsFunctionExpr<> y("y",2);
  variable fx = A.getCoeff(x);
  variable fy = A.getCoeff(y);
  auto j00 = fjac(fx).tr()*jac(G)*fjac(fx);
  auto j10 = fjac(fy).tr()*jac(G)*fjac(fx);
  auto j01 = fjac(fx).tr()*jac(G)*fjac(fy);
  auto j11 = fjac(fy).tr()*jac(G)*fjac(fy);

  auto detJ = j00*j11 - j01*j10;

  A.assemble(u*detJ);

  solVector = solver.solve(A.rhs());

  // for (index_t i = 0; i < solVector.size(); i++){
  //     if (solVector(i,0) > 0){
  //       gsInfo << "\n i = " << i << "is a problematic point...\n";
  //       solVector(i,0) = 1;
  //   } else {
  //       solVector(i,0) = 0;
  //   }
  // }

  gsMultiPatch<> dJ;
  u_sol.extract(dJ);

  variable out = A.getCoeff(dJ);

	gsInfo<<"Plotting " << name << " in Paraview...\n";
  ev.writeParaview( out   , G, name);
	// ev.options().setSwitch("plot.elements", true);

}

// void detJacConstraint::getDerivRhsFromPatch(index_t patch, gsSparseMatrix<> &xJac, gsSparseMatrix<> &yJac){
//     gsMultiPatch<> singlePatch(mp->patch(patch));
//
//     gsMultiBasis<> dJbas(m_detJacBasis.basis(patch));
//
//     // Prepare assembler
//     gsExprAssembler<> A(1,1);
//     gsMultiBasis<> dbasis(singlePatch);
//     A.setIntegrationElements(dJbas);
//     gsExprEvaluator<> ev(A);
//
//     // Define types
//     typedef gsExprAssembler<>::geometryMap geometryMap;
//     typedef gsExprAssembler<>::variable    variable;
//     typedef gsExprAssembler<>::space       space;
//     typedef gsExprAssembler<>::solution    solution;
//
//     geometryMap G = A.getMap(singlePatch);
//
//     // Use simple functions to generate standard basis vectors [1,0] and [0,1]
//     gsFunctionExpr<> x("x",2);
//     gsFunctionExpr<> y("y",2);
//     variable fx = ev.getVariable(x);
//     variable fy = ev.getVariable(y);
//
//     // Setup and assemble the two matricies
//     space u = A.getSpace(dbasis);
//     space v = A.getTestSpace(u,dJbas);
//
//     A.initMatrix();
//     auto j00 = grad(fx)*jac(G)*grad(fx).tr();
//     auto j10 = grad(fy)*jac(G)*grad(fx).tr();
//     auto j01 = grad(fx)*jac(G)*grad(fy).tr();
//     auto j11 = grad(fy)*jac(G)*grad(fy).tr();
//
//     auto dxiR = grad(u)*grad(fx).tr();
//     auto detaR = grad(u)*grad(fy).tr();
//
//     auto v1 = dxiR*j11;
//     auto v2 = -detaR*j10;
//
//     auto z0gradu = v1 + v2;
//     A.assemble(v*z0gradu.tr());
//
//     xJac = A.matrix();
//
//     A.initSystem();
//
//     auto v3 = -dxiR*j01;
//     auto v4 = detaR*j00;
//
//     auto z1gradu = v3 + v4;
//     A.assemble(v*z1gradu.tr());
//
//     yJac = A.matrix();
// }
