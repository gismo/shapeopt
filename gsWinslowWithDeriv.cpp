#include <gismo.h>
#include "gsWinslowWithDeriv.h"
using namespace gismo;

gsWinslowWithDeriv::gsWinslowWithDeriv(memory::shared_ptr<gsMultiPatch<>> mpin, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps):
    gsWinslow(mpin,use_dJC,useTensorStructureforDJC,checkForInf,checkForInf_eps)
    //, m_aff(this,false)
{
    // m_aff.computeMap();
}

gsWinslowWithDeriv::gsWinslowWithDeriv(memory::shared_ptr<gsMultiPatch<>> mpin, std::vector< gsDofMapper > mappers, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps):
    gsWinslow(mpin,mappers,use_dJC,useTensorStructureforDJC,checkForInf,checkForInf_eps)
    //, m_aff(this,false)
{
    // m_aff.computeMap();
}

// Returns false if the update failed
bool gsWinslowWithDeriv::update(gsVector<> x)
{
    gsVector<> tagged = getTagged();
    gsVector<> flat = getFlat();
    real_t diff = (x-tagged).norm();
    // gsMatrix<> mat(n_tagged,2);
    // mat << x,tagged;
    // gsInfo << mat <<"\n\n";

    gsInfo << "Difference in tagged and x: " << diff << "\n";

    if (diff == 0){
        return true;
    }

    // m_aff.update(x);

    updateTagged(x);
    m_aff.update();

    // Check for negative determinant before updating
    if (!checkForNegativeDetJ()){
        gsInfo << "Negative jac found\n";
        updateFlat(flat); // go back to original design
        return false;
    }

    // gsInfo << "NO NEG DET JAC\n";

    // gsWriteParaview(*m_mp,"geom");

    bool status = update();

    m_aff.reset();
    // m_aff.computeMap();

    return status;

}

// Note that this can be implemented more efficient, if we don't use m_affOPM
gsMatrix<> gsWinslowWithDeriv::jacobUpdate(gsVector<> x)
{
    // FIXIT: Check if we already updated !
    update(x);

    gsInfo << "Norm of grad = " << gradObj().norm() << "\n";

    gsMatrix<> Hcx;
    gsMatrix<> Hcc = hessObj(Hcx);

    Eigen::FullPivLU<Eigen::Matrix<real_t,Dynamic,Dynamic>> solver;
    solver.compute(Hcc);

    gsMatrix<> jac = solver.solve(-Hcx);

    return jac;
}

bool gsWinslowWithDeriv::checkForNegativeDetJ(){

    // gsInfo << "norm of flat : " << getFlat().norm() << "\n";
    // gsInfo << "norm of tagged : " << getTagged().norm() << "\n";
    // gsInfo << "quA, quB = " << m_quA << ", " << m_quB << "\n";
    // gsInfo << "MIN JAC G = " << minDJ << "\n";
    real_t minDJ = minDetJInGaussPts();

    if (minDJ <= m_checkForInf_eps) // Check for too small value of detJ at gauss points
    {
        return false;
    }

    return true;
}

real_t gsWinslowWithDeriv::minDetJInGaussPts(index_t incPts){

    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB + incPts);
    ev.options().setReal("quA",m_quA);

    typedef gsExprAssembler<>::geometryMap geometryMap;

    geometryMap G = A.getMap(*m_mp);
	gsDebugVar(ev.max((jac(G).det())));

    return ev.min(jac(G).det());
}
