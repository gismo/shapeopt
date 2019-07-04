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
    // by setting the degree to (2p-1) for 2D, and (3p-1) for 3D
    int p = m_detJacBasis.maxCwiseDegree();
    m_detJacBasis.setDegree(m_mp->targetDim()*p-1);
    // and reducing the continuity
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

    // FIXIT: Clean this up
    // Save mapper
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(m_detJacBasis);
    gsExprEvaluator<> ev(A);

    // Define types
    typedef gsExprAssembler<>::space space;
    space v = A.getSpace(dbasis);

    A.initSystem();
    m_space_mapper = v.mapper();

    // gsInfo << "n_controlpoints = " << n_controlpoints << "\n \n";
    // gsInfo << "n_constraints = " << n_controlpoints << "\n \n";
}

void gsDetJacConstraint::evalCon_into(gsAsVector<real_t> & result)
{
    result = evalCon();
}

gsVector<> gsDetJacConstraint::evalCon(){
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    index_t start = 0;
    gsVector<> result(m_size);
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

        // gsJacDetField<real_t> jacDetField(singlePatch);
        // variable detJ = ev.getVariable(jacDetField);
        // A.assemble(u*detJ);

        gsMultiPatch<> singlePatch(m_mp->patch(i));
        geometryMap G = A.getMap(singlePatch);
        A.assemble(u*jac(G).det());

        solVector = m_solversMassMatrix[i].solve(A.rhs());

        // Save result in result vector
        index_t len   = m_detJacBasis.size(i);
        result.segment(start,len) = solVector;

        start += m_detJacBasis.size(i);
    }

    return result;

}

gsSparseMatrix<> gsDetJacConstraint::getDerivRhsFromPatch(index_t patch){
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

    // gsInfo << "\nsize (" << A.matrix().rows() << ", " << A.matrix().cols() << ")\n";
    // index_t r = A.matrix().rows();
    // index_t c = A.matrix().cols();
    // xJac = A.matrix().block(0,0,r,c/2);
    // yJac = A.matrix().block(0,c/2,r,c/2);

    return A.matrix();


}

gsMatrix<> gsDetJacConstraint::getJacobianFromPatch(index_t patch){
    GISMO_ASSERT(m_areSolversSetup[patch],"Solver is not setup before calling detJacConstraint::getJacobianFromPatch");

    gsSparseMatrix<> Rhs = getDerivRhsFromPatch(patch);

    return m_solversMassMatrix[patch].solve(Rhs);
}

gsIpOptSparseMatrix gsDetJacConstraint::getJacobian(){
    // For each patch generate Jacobian with respect to all coordinates of geometry
    index_t dim = m_mp->targetDim();

    // Store a vector of pointers to the object, to avoid calling the constructor
    std::vector< memory::unique_ptr<gsIpOptSparseMatrix> > vMat(3);
    for(index_t i = 0; i < m_mp->nBoxes(); i++){ // For each patch
        // get xJac and yJac (and zJac in 3D) for patch i
        gsMatrix<> jac = getJacobianFromPatch(i);
        index_t r = jac.rows();
        index_t c = jac.cols()/dim;

        for (index_t d = 0; d < dim; d++)
        {
            gsMatrix<> mat = jac.block(0,d*c,r,c); // FIXIT: allow gsIpOptSparseMatrix to load this matrix directly by using const matrix as input
            if (i == 0){
                vMat[d].reset(new gsIpOptSparseMatrix(mat,-1));    // -1 indicates that IpOptSparseMatrix should generate dense matrix
            } else {
                vMat[d]->concatenate(gsIpOptSparseMatrix(mat,-1),"diag");
            }
        }
    }

    // Store full jacobian in first element of vMat
    for (index_t d = 1; d < dim; d++)
    {
        vMat[0]->concatenate(*vMat[d],"row");
    }

    return *vMat[0];

}

gsVector<> gsDetJacConstraint::getUpperBounds(){
    gsVector<> out;
    out.setConstant(n_constraints, 1e19);
    return out;
}

gsVector<> gsDetJacConstraint::getLowerBounds(){
    gsVector<> out;
    out.setConstant(n_constraints, m_eps);
    return out;
}

gsMultiPatch<> gsDetJacConstraint::getDetJ(){
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
    return dJ;
}

void gsDetJacConstraint::plotDetJ(std::string name){
    gsMultiPatch<> dJ = getDetJ();

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

    variable out = A.getCoeff(dJ);

    gsInfo<<"Plotting " << name << " in Paraview...\n";
    ev.writeParaview( out   , G, name);
    ev.options().setSwitch("plot.elements", true);

}

gsSparseMatrix<> gsDetJacConstraint::getMassMatrix(index_t i)
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    index_t start = 0;
    gsVector<> result(m_size);
    gsExprAssembler<> A(1,1);

    gsMultiBasis<> dbasis(m_detJacBasis.basis(i));
    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    space u = A.getSpace(dbasis);

    A.initSystem();
    A.assemble(u*u.tr());
    return A.matrix();
}

index_t gsDetJacConstraint::getSignOfPatch(index_t patch){
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsMultiBasis<> dbasis(m_detJacBasis.basis(patch));
    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    space u = A.getSpace(dbasis);

    A.initSystem();
    if (! m_areSolversSetup[patch]){
        A.assemble(u*u.tr());
        m_solversMassMatrix[patch].compute(A.matrix());
        m_areSolversSetup[patch] = true;
    }

    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    gsMultiPatch<> singlePatch(m_mp->patch(patch));
    geometryMap G = A.getMap(singlePatch);
    A.assemble(u*jac(G).det());

    solVector = m_solversMassMatrix[patch].solve(A.rhs());

    real_t avg = solVector.sum()/solVector.size();

    return (avg > 0) - (avg < 0); // returning sign of avg

}
