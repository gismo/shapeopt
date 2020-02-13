#include <gismo.h>
#include "gsConstraint.h"
#include "gsDetJacConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

gsDetJacConstraint::gsDetJacConstraint(memory::shared_ptr<gsMultiPatch<>> mpin, bool useTPSolver):
    gsConstraint(mpin),
    m_solversMassMatrix(mpin->nBoxes()),
    m_areSolversSetup(mpin->nBoxes()),
    m_useTPSolver(useTPSolver),
    m_solversTensor(mpin->nBoxes())
{
    // GISMO_ERROR("STOP");
    setup(); // Method to setup m_detJacBasis and m_space_mapper

    // gsInfo << "n_controlpoints = " << n_controlpoints << "\n \n";
    // gsInfo << "n_constraints = " << n_controlpoints << "\n \n";
}

gsVector<> gsDetJacConstraint::evalCon()
{
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
            if (m_useTPSolver){
                if (m_mp->domainDim() == 2){
                    m_solversTensor[i] = gsPatchPreconditionersCreator<>::massMatrixInvOp(static_cast<gsTensorBSplineBasis<2>&> (m_detJacBasis.basis(i)));
                } else { // If dimension is not 2D assume 3D.
                    m_solversTensor[i] = gsPatchPreconditionersCreator<>::massMatrixInvOp(static_cast<gsTensorBSplineBasis<3>&> (m_detJacBasis.basis(i)));
                }
            } else {
                A.assemble(u*u.tr());
                m_solversMassMatrix[i].compute(A.matrix());
            }
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

        if (m_useTPSolver){
            m_solversTensor[i]->apply(A.rhs(), solVector);
        } else {
            solVector = m_solversMassMatrix[i].solve(A.rhs());
        }

        // Save result in result vector
        index_t len   = m_detJacBasis.size(i);
        result.segment(start,len) = solVector;

        start += m_detJacBasis.size(i);
    }

    return result;

}

gsSparseMatrix<> gsDetJacConstraint::getDerivRhsFromPatch(index_t patch)
{
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

gsMatrix<> gsDetJacConstraint::getJacobianFromPatch(index_t patch)
{
    GISMO_ASSERT(m_areSolversSetup[patch],"Solver is not setup before calling detJacConstraint::getJacobianFromPatch");

    gsSparseMatrix<> Rhs = getDerivRhsFromPatch(patch);

    if (m_useTPSolver){
        gsMatrix<> out;
        m_solversTensor[patch]->apply(Rhs,out);
        return out;
    } else {
        return m_solversMassMatrix[patch].solve(Rhs);
    }
}

gsIpOptSparseMatrix gsDetJacConstraint::getJacobian()
{
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

gsVector<> gsDetJacConstraint::getUpperBounds()
{
    gsVector<> out;
    out.setConstant(n_constraints, 1e19);
    return out;
}

gsVector<> gsDetJacConstraint::getLowerBounds()
{
    gsVector<> out;
    out.setConstant(n_constraints, m_eps);
    return out;
}

gsMultiPatch<> gsDetJacConstraint::getDetJ()
{

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    // Elements used for numerical integration
    A.setIntegrationElements(m_detJacBasis);
    gsExprEvaluator<> ev(A);

    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);

    space u = A.getSpace(m_detJacBasis);

    geometryMap G = A.getMap(*m_mp);

    A.initSystem();

    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    A.assemble(u*jac(G).det());

    if (m_useTPSolver && m_mp->nBoxes() == 1)
        m_solversTensor[0]->apply(A.rhs(), solVector);
    else
    {
        A.assemble(u*u.tr());
        gsSparseSolver<>::LU solver;
        solver.compute(A.matrix());
        solVector = solver.solve(A.rhs());
        // solVector = evalCon();
    }

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

void gsDetJacConstraint::plotDetJ(std::string name)
{
    gsMultiPatch<> dJ = getDetJ();

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);

    // Elements used for numerical integration
    A.setIntegrationElements(m_detJacBasis);
    gsExprEvaluator<> ev(A);

    space u = A.getSpace(m_detJacBasis);

    geometryMap G = A.getMap(*m_mp);

    variable out = A.getCoeff(dJ);

    gsInfo<<"Plotting " << name << " in Paraview...\n";
    // gsWriteParaview(dJ,name,100000,true);
    ev.options().setSwitch("plot.elements", true);
    ev.options().setInt("plot.npts", 50000);
    ev.writeParaview( out   , G, name);

}

gsMultiPatch<> gsDetJacConstraint::getDetJFromCoef()
{
    GISMO_ASSERT(m_mp->targetDim() == 2, "getDetJSurface only works for 2D");

    // gsMultiPatch in which to store the surface
    gsMultiPatch<> detJ;

    gsVector<> d = evalCon();
    index_t start = 0;

    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        gsMultiPatch<> singlePatch(m_mp->patch(p));
        gsMultiBasis<> multibas(m_detJacBasis.basis(p));

        // Prepare coefficient matrix
        gsMatrix<> greville = m_detJacBasis.basis(p).anchors();
        gsMatrix<> coefs;
        coefs.setOnes(greville.cols(), 1);

        index_t len   = m_detJacBasis.size(p);
        coefs.col(0) = d.segment(start,len);
        start += len;

	    gsTHBSpline<2> geom(m_detJacBasis.basis(p), coefs);
        detJ.addPatch(geom);

    }

    return detJ;



}

gsMultiPatch<> gsDetJacConstraint::getDetJSurface(bool zero)
{
    GISMO_ASSERT(m_mp->targetDim() == 2, "getDetJSurface only works for 2D");

    // gsMultiPatch in which to store the surface
    gsMultiPatch<> detJSurface;

    gsVector<> d = evalCon();
    index_t start = 0;

    real_t max = 1;
    for (index_t i = 0; i < d.size(); i++)
    {
        if (d[i] > max)
            d[i] = max;

        if (zero)
            d[i] = 0;
    }

    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        gsMultiPatch<> singlePatch(m_mp->patch(p));
        gsMultiBasis<> multibas(m_detJacBasis.basis(p));

        // Prepare coefficient matrix
        gsMatrix<> greville = m_detJacBasis.basis(p).anchors();
        gsMatrix<> coefs;
        coefs.setOnes(greville.cols(), 3);

        index_t len   = m_detJacBasis.size(p);
        coefs.col(2) = d.segment(start,len);
        start += len;

        // Project G into m_detJacBasis's space
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;

        gsExprAssembler<> A(1,1);

        // Elements used for numerical integration
        A.setIntegrationElements(multibas);
        gsExprEvaluator<> ev(A);

        A.options().setInt("quB",m_quB);
        A.options().setReal("quA",m_quA);

        space u = A.getSpace(m_detJacBasis.basis(p));

        gsMatrix<> solVector;

        A.initSystem();

        for (index_t dim = 0; dim < m_mp->targetDim(); dim++)
        {
            gsMultiPatch<> x = singlePatch.coord(dim);
            variable xvar = ev.getVariable(x);

            A.initSystem();
            A.assemble(u*xvar);

            gsMatrix<> rhs = A.rhs();

            if (m_useTPSolver){
                m_solversTensor[p]->apply(rhs, solVector);
            } else {
                solVector = m_solversMassMatrix[p].solve(rhs);
            }

            coefs.col(dim) = solVector;
        }

        // gsGeometry<>::uPtr gg = gsNurbsCreator<>::BSplineSquare();
        // gsHBSplineBasis<2> hb(gg->basis());
        // gsMatrix<> mm(hb.size(),1); mm.setZero();
        // gsHBSpline<2> geom(hb,mm);

	    gsTHBSpline<2> geom(m_detJacBasis.basis(p), coefs);
        detJSurface.addPatch(geom);

    }

    return detJSurface;



}

real_t gsDetJacConstraint::refineDetJSurfaceUntilPositive(index_t nRefSteps, gsMultiPatch<> & dJ)
{

    // Does not work as long as solver is Tensor product
    GISMO_ASSERT(!m_useTPSolver,"Cannot use adaptive refinement when Tensor Product solver is enabled.\n");

    // Create gsMultiBasis to hold the adaptive refinable basis, and fill it with H-spline bases


    // Create copy of m_detJacBasis.
    gsMultiBasis<> basis(m_detJacBasis);

    // Clear m_detJacBasis to delete all bases inside,
    //      to make room for the gsHBSplineBasis version that allow adaptivity
    m_detJacBasis.clear();
    gsInfo << "Size of m_detJacBasis : " << m_detJacBasis.size() << "\n";

    // Fill it with H-spline bases made from those in basis
    for (index_t pn = 0; pn < basis.nPieces(); pn ++){
        gsTensorBSplineBasis<2,real_t> bas = dynamic_cast<gsTensorBSplineBasis<2,real_t>&>(basis.basis(pn));
        gsTHBSplineBasis<2> HBbas(bas);
        m_detJacBasis.addBasis(HBbas.clone());
    }

    real_t minD;

    for (index_t i = 0; i < nRefSteps; i++)
    {
        // Reset the class with the new basis
        reset();

        minD = evalCon().minCoeff();

        // If all coefficients, i.e. the minimal one, is positive stop the loop
        // if (minD > 0)
        //     return minD;

        // Mark elements where a basis function with negative coefficient has support
        std::vector<bool> elMarked;
        markElements(elMarked,0);
        refineElements(elMarked,m_detJacBasis);
        refineElements(elMarked,dJ);

    }

    testSplineDetJ();
    setup();
    return minD;
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

void gsDetJacConstraint::plotActiveConstraints(std::vector<bool> & elMarked, std::string name, real_t tol1)
{
    // Create gsMultiPatch to save sum of basis functions which supporting elements should be marked
    gsMultiPatch<> basMarking;
    for(int p = 0; p < m_mp->nBoxes(); p++)
    {
        gsMatrix<> coefs;
        coefs.setZero(m_detJacBasis.size(p),1);
        basMarking.addPatch( m_detJacBasis.basis(p).makeGeometry( coefs ) );
    }

    gsVector<> d = evalCon();

    // Loop through d-vector and mark
    index_t start = 0;
    for(int p = 0; p < m_mp->nBoxes(); p++){
        // Save result in result vector
        index_t len = m_detJacBasis.size(p);
        for(int i = 0; i < len; i++){
            if (d[start + i] - m_eps < tol1 ){
                // gsInfo << "d[" << start + i << "] = " << d[start + i] << " and m_eps = " << m_eps << " and tol1 = " << tol1 << "\n";
                basMarking.patch(p).coef(i,0) = 1;
            }
        }
        start += m_detJacBasis.size(p);
    }

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

    variable out = A.getCoeff(basMarking);

    gsInfo<<"Plotting " << name << " in Paraview...\n";
    ev.writeParaview( out   , G, name);
    ev.options().setSwitch("plot.elements", true);

}

void gsDetJacConstraint::setup()
{
    gsMultiBasis<> bas(*m_mp);
    m_detJacBasis = bas;

    m_Hspline_flag = false;

    // Prepare basis for detJac
    // by setting the degree to (2p-1) for 2D, and (3p-1) for 3D
    int p = m_detJacBasis.maxCwiseDegree();
    m_detJacBasis.setDegree(m_mp->targetDim()*p-1);
    // and reducing the continuity
    m_detJacBasis.reduceContinuity(1);
    // gsInfo << "\n..detJacBasis degree is: " << 2*p-1 << "\n";

    reset();

}

void gsDetJacConstraint::reset(){
    for (index_t i = 0; i < m_mp->nBoxes(); i++){
        m_areSolversSetup[i] = false;
    }

    m_size = m_detJacBasis.size();
    // gsInfo << "m_size : " << m_size << "\n";
    // gsInfo << "m_eps  : " << m_eps << "\n";

    // Count the total number of controlpoints
    n_controlpoints = 0;
    n_constraints = 0;
    for(int i = 0; i < m_mp->nBoxes(); i++){
        n_controlpoints += m_mp->patch(i).coefsSize();
        n_constraints += m_detJacBasis.size(i);
        // gsInfo << "Patch " << i << " has " << mp->patch(i).coefsSize() << " cps... \n";
        // gsInfo << "Patch " << i << " gets " << m_detJacBasis.size(i) << " constraints... \n";
    }
    // gsInfo << "n_constraints = " << n_constraints << "\n";

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
}

gsVector<> gsDetJacConstraint::massMatrixSolve(gsVector<> rhs) const
{
    index_t start = 0;
    gsVector<> result(m_size);
    for(int i = 0; i < m_mp->nBoxes(); i++){
        if (! m_areSolversSetup[i]){
            GISMO_ERROR("solvers are not setup before massMatrixSolve in gsDetJacConstraint\n");
        }

        gsMatrix<> solVector;
        index_t len   = m_detJacBasis.size(i);

        if (m_useTPSolver){
            m_solversTensor[i]->apply(rhs.segment(start,len), solVector);
        } else {
            solVector = m_solversMassMatrix[i].solve(rhs.segment(start,len));
        }

        // Save result in result vector
        result.segment(start,len) = solVector;

        start += m_detJacBasis.size(i);
    }
    return result;
}

gsSparseMatrix<> gsDetJacConstraint::hessD(index_t patch, index_t i) const
{
    // A hack to get a basis function
    gsMultiPatch<> bas;
    gsMatrix<> coefs;
    for(index_t p = 0; p < m_mp->nBoxes(); p++){
        coefs.setZero(m_detJacBasis.size(p),1);
        if (p == patch)
            coefs(i,0) = 1;

        bas.addPatch( m_detJacBasis.basis(p).makeGeometry( coefs ));

    }

    // Prepare assembler
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(m_detJacBasis);
    gsExprEvaluator<> ev(A);

    // Define types
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    // Setup and assemble the two matricies
    space u = A.getSpace(dbasis,m_mp->geoDim());

    variable b = ev.getVariable(bas);

    A.initSystem();

    auto ddetJdc =  (jac(u)%jac(G).inv().tr())*jac(G).det();
    auto dJinvTdc = - matrix_by_space_tr(jac(G).inv(),jac(u))*jac(G).inv().tr();

    A.assemble(
        ddetJdc*(b.val()*matrix_by_space(jac(G).inv(),jac(u)).trace()).tr()
        +
        (jac(u) % (b.val()*jac(G).det().val()*dJinvTdc))
    );


    return A.matrix();

}

index_t gsDetJacConstraint::sizeOfBasis(index_t k) const{
    return m_detJacBasis.size(k);
}

real_t gsDetJacConstraint::provePositivityOfDetJ_TP(index_t & neededRefSteps, index_t maxRefSteps)
{

    // Save old basis

    real_t minD;

    for (index_t i = 0; i < maxRefSteps; i++)
    {
        neededRefSteps = i;

        minD = evalCon().minCoeff();

        // If all coefficients, i.e. the minimal one, is positive stop the loop
        if (minD > 0)
        {
            setup();
            return minD;
        }

        // refine and setup the class with the new basis
        m_detJacBasis.uniformRefine();
        reset();
    }

    gsInfo << "DetJ could not be proven positive in " << maxRefSteps << "iterations. Smallest coef was still " << minD << "\n";
    setup();
    return minD;
}

real_t gsDetJacConstraint::refineUntilPositive(index_t maxRefSteps, real_t tol)
{

    // Does not work as long as solver is Tensor product
    GISMO_ASSERT(!m_useTPSolver,"Cannot use adaptive refinement when Tensor Product solver is enabled.\n");

    // Create gsMultiBasis to hold the adaptive refinable basis, and fill it with H-spline bases

    if (not m_Hspline_flag)
    {
        // Create copy of m_detJacBasis.
        gsMultiBasis<> basis(m_detJacBasis);

        // Clear m_detJacBasis to delete all bases inside,
        //      to make room for the gsHBSplineBasis version that allow adaptivity
        m_detJacBasis.clear();
        // gsInfo << "Size of m_detJacBasis : " << m_detJacBasis.size() << "\n";

        // Fill it with H-spline bases made from those in basis
        for (index_t pn = 0; pn < basis.nPieces(); pn ++){
            gsTensorBSplineBasis<2,real_t> bas = dynamic_cast<gsTensorBSplineBasis<2,real_t>&>(basis.basis(pn));
            gsTHBSplineBasis<2> HBbas(bas);
            m_detJacBasis.addBasis(HBbas.clone());
        }

        m_Hspline_flag = true;
    }


    // gsInfo << "Size of m_detJacBasis : " << m_detJacBasis.size() << "\n";

    real_t minD;

    for (index_t i = 0; i < maxRefSteps; i++)
    {
        // Reset the class with the new basis
        reset();

        minD = evalCon().minCoeff();

        // If all coefficients, i.e. the minimal one, is positive stop the loop
        if (minD > tol)
            return minD;

        // Mark elements where a basis function with negative coefficient has support
        std::vector<bool> elMarked;
        markElements(elMarked,tol);
        refineElements(elMarked,m_detJacBasis);

    }

    reset();
    minD = evalCon().minCoeff();
    gsInfo << "DetJ could not be proven positive in " << maxRefSteps << "iterations. Smallest coef was still " << minD << "\n";
    return minD;
}

// FIXIT: DEBUG
real_t gsDetJacConstraint::provePositivityOfDetJ(index_t maxRefSteps)
{

    // Does not work as long as solver is Tensor product
    GISMO_ASSERT(!m_useTPSolver,"Cannot use adaptive refinement when Tensor Product solver is enabled.\n");

    // Save old basis
    gsMultiBasis<> tmpBasis(m_detJacBasis);

    // Refine until positive detJ coefs
    real_t minD = refineUntilPositive(maxRefSteps);

    // Go back to original basis
    m_detJacBasis = tmpBasis;
    return minD;
}

// FIXIT: DEBUG
void gsDetJacConstraint::refineElements(std::vector<bool> & elMarked, gsMultiBasis<> & basis)
{

    const int dim = basis.dim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    int numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<> refBoxes;


    for (unsigned pn=0; pn < basis.nBases(); ++pn )// for all patches
    {
        // Get number of elements to be refined on this patch
        const int numEl = basis[pn].numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  std::bind2nd(std::equal_to<bool>(), true) );
        gsInfo << "Refine " << numMarked << "elements\n" << std::flush;
        poffset += numEl;
        refBoxes.resize(dim, 2*numMarked);
        //gsDebugVar(numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<>::domainIter domIt = basis.basis(pn).makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount++ ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                refBoxes.col(2*numMarked  ) =
                        refBoxes.col(2*numMarked+1) = domIt->centerPoint();

                // Advance marked cells counter
                numMarked++;
            }
        }
        // Refine all of the found refBoxes in this patch
        (dynamic_cast< gsTHBSplineBasis<2> &>(basis.basis(pn))).refine( refBoxes );
        // Check if basis.basis(pn) is refined..

    }
    // gsInfo << "Size of basis = " << basis.size() << "\n";
}

// FIXIT: Implement refExtension functionally.. ask Angelos
void gsDetJacConstraint::refineElements(const std::vector<bool> & elMarked, gsMultiPatch<> & mp)
{
    gsMultiBasis<> basis(mp);
    const int dim = basis.dim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    int numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<> refBoxes;

    for (unsigned pn=0; pn < basis.nBases(); ++pn )// for all patches
    {
        // Get number of elements to be refined on this patch
        const int numEl = basis[pn].numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  std::bind2nd(std::equal_to<bool>(), true) );
        poffset += numEl;
        refBoxes.resize(dim, 2*numMarked);
        //gsDebugVar(numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<>::domainIter domIt = basis.basis(pn).makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount++ ] ) // refine this element ?
            {
                // Construct degenerate box by setting both
                // corners equal to the center
                refBoxes.col(2*numMarked  ) =
                        refBoxes.col(2*numMarked+1) = domIt->centerPoint();

                // Advance marked cells counter
                numMarked++;
            }
        }
        // Refine all of the found refBoxes in this patch
        // FIXIT: make dimension independet.. Find way to replace 2 with dim.
        std::vector<gsSortedVector<unsigned> > OX = dynamic_cast<gsHTensorBasis<2,real_t>&> (mp.patch(pn).basis()).getXmatrix();
        mp.patch(pn).basis().refine( refBoxes );
        gsSparseMatrix<> transf;
        dynamic_cast<gsHTensorBasis<2,real_t>&> (mp.patch(pn).basis()).transfer(OX, transf);
        //gsDebug<<"tranf orig:\n"<<transf<<std::endl;
        mp.patch(pn).coefs() = transf*mp.patch(pn).coefs();

    }

    mp.repairInterfaces();
}

void gsDetJacConstraint::markElements(std::vector<bool> & elMarked, gsVector<> lambda, real_t tol)
{
    // Create gsMultiPatch to save sum of basis functions which supporting elements should be marked
    gsMultiPatch<> basMarking;
    for(int p = 0; p < m_mp->nBoxes(); p++)
    {
        gsMatrix<> coefs;
        coefs.setZero(m_detJacBasis.size(p),1);
        basMarking.addPatch( m_detJacBasis.basis(p).makeGeometry( coefs ) );
    }

    // Loop through lambda and mark
    index_t start = 0;
    for(int p = 0; p < m_mp->nBoxes(); p++){
        // Save result in result vector
        index_t len = m_detJacBasis.size(p);
        for(int i = 0; i < len; i++){
            if (lambda[start + i] > tol ){
                gsInfo << "lambda[" << start + i << "] = " << lambda[start + i] << " and m_eps = " << m_eps << " and tol = " << tol << "\n" << std::flush;
                basMarking.patch(p).coef(i,0) = 1;
            }
        }
        start += m_detJacBasis.size(p);
    }

    // Integrate element wise to get support marked
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(m_detJacBasis);
    gsExprEvaluator<>::variable markfun = ev.getVariable(basMarking);
    // Get the element-wise norms.
    ev.integralElWise(markfun);
    const std::vector<real_t> & eltErrs  = ev.elementwise();

    // Mark with true false
    elMarked.resize( eltErrs.size() );
    // Now just check for each element, whether the local error
    // is above the computed threshold or not, and mark accordingly.

    typename std::vector<real_t>::const_iterator err = eltErrs.begin();
    for(std::vector<bool>::iterator i = elMarked.begin(); i!=  elMarked.end(); ++i, ++err)
        *i = ( *err > 0 );
}

void gsDetJacConstraint::markElements(std::vector<bool> & elMarked, real_t tol1 , real_t tol2 )
{
    // Create gsMultiPatch to save sum of basis functions which supporting elements should be marked
    gsMultiPatch<> basMarking;
    for(int p = 0; p < m_mp->nBoxes(); p++)
    {
        gsMatrix<> coefs;
        coefs.setZero(m_detJacBasis.size(p),1);
        basMarking.addPatch( m_detJacBasis.basis(p).makeGeometry( coefs ) );
    }

    gsVector<> d = evalCon();

    // Loop through d-vector and mark
    index_t start = 0;
    for(int p = 0; p < m_mp->nBoxes(); p++){
        // Save result in result vector
        index_t len = m_detJacBasis.size(p);
        for(int i = 0; i < len; i++){
            if (d[start + i] < tol1 ){
                gsInfo << "d[" << start + i << "] = " << d[start + i] << " and m_eps = " << m_eps << " and tol1 = " << tol1 << "\n";
                basMarking.patch(p).coef(i,0) = 1;
            }
        }
        start += len;
    }

    // Integrate element wise to get support marked
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(m_detJacBasis);
    gsExprEvaluator<>::variable markfun = ev.getVariable(basMarking);
    // Get the element-wise norms.
    ev.integralElWise(markfun);
    const std::vector<real_t> & eltErrs  = ev.elementwise();

    // Mark with true false
    elMarked.resize( eltErrs.size() );
    // Now just check for each element, whether the local error
    // is above the computed threshold or not, and mark accordingly.

    typename std::vector<real_t>::const_iterator err = eltErrs.begin();
    for(std::vector<bool>::iterator i = elMarked.begin(); i!=  elMarked.end(); ++i, ++err)
        *i = ( *err > tol2 );

}

void gsDetJacConstraint::testSplineDetJ()
{
    gsMultiPatch<> dJ = getDetJ();
    testSplineDetJ(dJ);
}

void gsDetJacConstraint::testSplineDetJ(gsMultiPatch<> & dJ)
{

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);

    // Elements used for numerical integration
    A.setIntegrationElements(m_detJacBasis);
    gsExprEvaluator<> ev(A);

    space u = A.getSpace(m_detJacBasis);

    geometryMap G = A.getMap(*m_mp);

    variable dj = ev.getVariable(dJ);

    real_t diff = ev.integral(dj.val() - jac(G).det());
    gsInfo << "diff = " << diff << "\n";
}
