#include <gismo.h>
#include "gsConstraint.h"
#include "gsNewDetJacConstraint.h"
#include "gsIpOptSparseMatrix.h"
using namespace gismo;

gsNewDetJacConstraint::gsNewDetJacConstraint(memory::shared_ptr<gsMultiPatch<>> mpin):
    gsConstraint(mpin),
    m_mappers(mpin->targetDim())
{
    setup();
}

// Implement
gsVector<> gsNewDetJacConstraint::evalCon()
{
    // ... put code here ...
    GISMO_NO_IMPLEMENTATION;
}


// Implement
gsIpOptSparseMatrix gsNewDetJacConstraint::getJacobian()
{
    // ... put code here ...
    GISMO_NO_IMPLEMENTATION;
}

gsVector<> gsNewDetJacConstraint::getUpperBounds()
{
    gsVector<> out;
    out.setConstant(n_constraints, 1e19);
    return out;
}

gsVector<> gsNewDetJacConstraint::getLowerBounds()
{
    gsVector<> out;
    out.setConstant(n_constraints, m_eps);
    return out;
}

void gsNewDetJacConstraint::setSpaceMapper(){
    // FIXIT: Clean this up
    // Save mapper
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    gsExprEvaluator<> ev(A);

    A.setIntegrationElements(dbasis);

    // Define types
    typedef gsExprAssembler<>::space space;
    space v = A.getSpace(dbasis);

    A.initSystem();
    m_space_mapper = v.mapper();
}

// Implement
void gsNewDetJacConstraint::setup()
{
    // Count the total number of controlpoints
    n_controlpoints = 0;
    for(int i = 0; i < m_mp->nBoxes(); i++){
        n_controlpoints += m_mp->patch(i).coefsSize();
    }

    setSpaceMapper();

    // Get basis for partial derivatives
    std::vector< gsMultiBasis<> > bases;
    for (index_t d = 0; d < m_mp->targetDim(); d++){
        gsMultiBasis<> bas(*m_mp);
        //FIXIT: check if space is correct
        gsInfo << "CHECK IF SPACE IS RIGHT\n";
        bas.degreeDecrease(1,d);    // Decrease degree by one in direction d
        bases.push_back(bas);
    }

    // Only works for 2D.
    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    // Define types
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    A.setIntegrationElements(bases[0]);

    // One point quadrature
    A.options().setInt("quB",1);
    A.options().setReal("quA",0);

    space u = A.getSpace(bases[0]);
    space v = A.getTestSpace(u,bases[1],1);

    A.initSystem();
    m_mappers[0] = u.mapper();
    m_mappers[1] = v.mapper();
    A.assemble(v*u.tr());

    n_constraints = A.matrix().nonZeros();
    gsInfo << "n_constraints = " << n_constraints << "\n";

    // Save in datastructure
    m_overlapping.setZero(m_mp->targetDim(),n_constraints);
    index_t iter = 0;

    for (int k=0; k<A.matrix().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A.matrix(),k); it; ++it)
    {
        m_overlapping(0,iter) = it.row();
        m_overlapping(1,iter) = it.col();
        iter++;
    }

    gsInfo << iter << "\n";
}

// FIXIT, implement for this method. OR BETTER implement generally in gsConstraint

// void gsNewDetJacConstraint::plotActiveConstraints(std::vector<bool> & elMarked, std::string name, real_t tol1)
// {
//     // Create gsMultiPatch to save sum of basis functions which supporting elements should be marked
//     gsMultiPatch<> basMarking;
//     for(int p = 0; p < m_mp->nBoxes(); p++)
//     {
//         gsMatrix<> coefs;
//         coefs.setZero(m_detJacBasis.size(p),1);
//         basMarking.addPatch( m_detJacBasis.basis(p).makeGeometry( coefs ) );
//     }
//
//     gsVector<> d = evalCon();
//
//     // Loop through d-vector and mark
//     index_t start = 0;
//     for(int p = 0; p < m_mp->nBoxes(); p++){
//         // Save result in result vector
//         index_t len = m_detJacBasis.size(p);
//         for(int i = 0; i < len; i++){
//             if (d[start + i] - m_eps < tol1 ){
//                 gsInfo << "d[" << start + i << "] = " << d[start + i] << " and m_eps = " << m_eps << " and tol1 = " << tol1 << "\n";
//                 basMarking.patch(p).coef(i,0) = 1;
//             }
//         }
//         start += m_detJacBasis.size(p);
//     }
//
//     typedef gsExprAssembler<>::geometryMap geometryMap;
//     typedef gsExprAssembler<>::variable    variable;
//     typedef gsExprAssembler<>::space       space;
//     typedef gsExprAssembler<>::solution    solution;
//
//     gsExprAssembler<> A(1,1);
//
//     // Elements used for numerical integration
//     A.setIntegrationElements(m_detJacBasis);
//     gsExprEvaluator<> ev(A);
//
//     space u = A.getSpace(m_detJacBasis);
//
//     geometryMap G = A.getMap(*m_mp);
//
//     variable out = A.getCoeff(basMarking);
//
//     gsInfo<<"Plotting " << name << " in Paraview...\n";
//     ev.writeParaview( out   , G, name);
//     ev.options().setSwitch("plot.elements", true);
//
// }

// FIXIT, implement for this method. OR BETTER implement generally in gsConstraint

// void gsNewDetJacConstraint::markElements(std::vector<bool> & elMarked, real_t tol1 , real_t tol2 )
// {
//     // Create gsMultiPatch to save sum of basis functions which supporting elements should be marked
//     gsMultiPatch<> basMarking;
//     for(int p = 0; p < m_mp->nBoxes(); p++)
//     {
//         gsMatrix<> coefs;
//         coefs.setZero(m_detJacBasis.size(p),1);
//         basMarking.addPatch( m_detJacBasis.basis(p).makeGeometry( coefs ) );
//     }
//
//     gsVector<> d = evalCon();
//
//     // Loop through d-vector and mark
//     index_t start = 0;
//     for(int p = 0; p < m_mp->nBoxes(); p++){
//         // Save result in result vector
//         index_t len = m_detJacBasis.size(p);
//         for(int i = 0; i < len; i++){
//             if (d[start + i] - m_eps < tol1 ){
//                 gsInfo << "d[" << start + i << "] = " << d[start + i] << " and m_eps = " << m_eps << " and tol1 = " << tol1 << "\n";
//                 basMarking.patch(p).coef(i,0) = 1;
//             }
//         }
//         start += m_detJacBasis.size(p);
//     }
//
//     // Integrate element wise to get support marked
//     gsExprEvaluator<> ev;
//     ev.setIntegrationElements(m_detJacBasis);
//     gsExprEvaluator<>::variable markfun = ev.getVariable(basMarking);
//     // Get the element-wise norms.
//     ev.integralElWise(markfun);
//     const std::vector<real_t> & eltErrs  = ev.elementwise();
//
//     // Mark with true false
//     elMarked.resize( eltErrs.size() );
//     // Now just check for each element, whether the local error
//     // is above the computed threshold or not, and mark accordingly.
//
//     typename std::vector<real_t>::const_iterator err = eltErrs.begin();
//     for(std::vector<bool>::iterator i = elMarked.begin(); i!=  elMarked.end(); ++i, ++err)
//         *i = ( *err > tol2 );
//
// }
