#include <gismo.h>
#include "gsDetJacConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

gsDetJacConstraint::gsDetJacConstraint(gsMultiPatch<>* mpin): m_mp(mpin), m_detJacBasis(*m_mp),
        m_solversMassMatrix(mpin->nBoxes()),m_areSolversSetup(mpin->nBoxes()){
    for (index_t i = 0; i < m_mp->nBoxes(); i++){
        m_areSolversSetup[i] = false;
    }
    // Prepare basis for detJac
    // by setting the degree to (2p-1)
    int p = m_detJacBasis.maxCwiseDegree();
    m_detJacBasis.setDegree(2*p-1);
    // by reducing the continuity
    m_detJacBasis.reduceContinuity(1);
    // gsInfo << "\n..detJacBasis degree is: " << 2*p-1 << "\n";
    m_size = m_detJacBasis.size();
    gsInfo << "m_size : " << m_size << "\n";

    // Count the total number of controlpoints
    n_controlpoints = 0;
    n_constraints = 0;
    for(int i = 0; i < m_mp->nBoxes(); i++){
        n_controlpoints += m_mp->patch(i).coefsSize();
        n_constraints += m_detJacBasis.size(i);
        // gsInfo << "Patch " << i << " has " << mp->patch(i).coefsSize() << " cps... \n";
        // gsInfo << "Patch " << i << " gets " << m_detJacBasis.size(i) << " constraints... \n";
    }

    // gsInfo << "n_controlpoints = " << n_controlpoints << "\n \n";
    // gsInfo << "n_constraints = " << n_controlpoints << "\n \n";
}

void gsDetJacConstraint::evalCon_into(gsVector<> &result){
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    index_t start = 0;
    for(int i = 0; i < m_mp->nBoxes(); i++){
        gsExprAssembler<> A(1,1);

        gsMultiBasis<> dbasis(m_detJacBasis.basis(i));
        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);

        space u = A.getSpace(dbasis);

        A.initSystem();
        if (! m_areSolversSetup[i]){
            A.assemble(u*u.tr());
            m_solversMassMatrix[i].compute(A.matrix());
            m_areSolversSetup[i] = true;
        }

        gsMatrix<> solVector;
        solution u_sol = A.getSolution(u, solVector);

        // gsJacDetField<real_t> jacDetField(mp->patch(i));
        // variable detJ = ev.getVariable(jacDetField);
        // A.assemble(u*detJ);

        geometryMap G = A.getMap(m_mp->patch(i));
        A.assemble(u*jac(G).det());

        solVector = m_solversMassMatrix[i].solve(A.rhs());

        // Save result in result vector
        index_t len   = m_detJacBasis.size(i);
        result.segment(start,len) = solVector;

        start += m_detJacBasis.size(i);
    }



}

gsVector<> gsDetJacConstraint::evalCon(){
    gsVector<> out(m_size);
    evalCon_into(out);
    return out;
}

void gsDetJacConstraint::getDerivRhsFromPatch(index_t patch, gsSparseMatrix<> &xJac, gsSparseMatrix<> &yJac){
    gsMultiPatch<> singlePatch(m_mp->patch(patch));

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
    space u = A.getSpace(dbasis,m_mp->geoDim());
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

void gsDetJacConstraint::getJacobianFromPatch(index_t patch, gsMatrix<> &xJac, gsMatrix<> &yJac){
    GISMO_ASSERT(m_areSolversSetup[patch],"Solver is not setup before calling detJacConstraint::getJacobianFromPatch");

    gsSparseMatrix<> xDrhs,yDrhs;
    getDerivRhsFromPatch(patch,xDrhs,yDrhs);

    xJac = m_solversMassMatrix[patch].solve(xDrhs);

    yJac = m_solversMassMatrix[patch].solve(yDrhs);
}

gsIpOptSparseMatrix gsDetJacConstraint::getJacobian(){
    // For each patch generate Jacobian with respect to x and y coordinates of geometry
    memory::unique_ptr<gsIpOptSparseMatrix> xMat,yMat;         // Store a pointer to the object, to avoid calling the constructor
    for(index_t i = 0; i < m_mp->nBoxes(); i++){
        // get xJac and yJac for patch i
        gsMatrix<> xJac,yJac;
        getJacobianFromPatch(i,xJac,yJac);
        if (i == 0){
            xMat.reset(new gsIpOptSparseMatrix(xJac,-1));    // -1 indicates that IpOptSparseMatrix should generate dense matrix
            yMat.reset(new gsIpOptSparseMatrix(yJac,-1));    // --||-- ...
        } else {
            xMat->concatenate(gsIpOptSparseMatrix(xJac,-1),"diag");
            yMat->concatenate(gsIpOptSparseMatrix(yJac,-1),"diag");
        }
    }
    xMat->concatenate(*yMat,"row");       // Store full jacobian in xMat

    return *xMat;

}

gsVector<> gsDetJacConstraint::getUpperBounds(){
    gsVector<> out;
    out.setConstant(n_constraints, -m_eps);
    return out;
}

gsVector<> gsDetJacConstraint::getLowerBounds(){
    gsVector<> out;
    out.setConstant(n_constraints, -1e19);
    return out;
}

void gsDetJacConstraint::plotDetJ(std::string name){
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

  geometryMap G = A.getMap(*m_mp);

  A.initSystem();
  A.assemble(u*u.tr());
  gsSparseSolver<>::CGDiagonal solver;
  solver.compute(A.matrix());

  gsMatrix<> solVector;
  solution u_sol = A.getSolution(u, solVector);

  A.assemble(u*jac(G).det());

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
