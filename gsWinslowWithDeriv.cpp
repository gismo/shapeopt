#include <gismo.h>
#include "gsWinslowWithDeriv.h"
using namespace gismo;

gsWinslowWithDeriv::gsWinslowWithDeriv(gsMultiPatch<>* mpin, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps):
    gsWinslow(mpin,use_dJC,useTensorStructureforDJC,checkForInf,checkForInf_eps),
    m_aff(this,false)
{
    // m_aff.computeMap();
}

gsWinslowWithDeriv::gsWinslowWithDeriv(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps):
    gsWinslow(mpin,mappers,use_dJC,useTensorStructureforDJC,checkForInf,checkForInf_eps),
    m_aff(this,false)
{
    // m_aff.computeMap();
}

// Returns false if the update failed
bool gsWinslowWithDeriv::update(gsVector<> x)
{
    gsVector<> tagged = getTagged();
    gsVector<> flat = getFlat();
    real_t diff = (x-tagged).norm();
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
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);

    typedef gsExprAssembler<>::geometryMap geometryMap;

    geometryMap G = A.getMap(*m_mp);

    real_t minDJ = ev.min(jac(G).det());

    // gsInfo << "norm of flat : " << getFlat().norm() << "\n";
    // gsInfo << "norm of tagged : " << getTagged().norm() << "\n";
    // gsInfo << "quA, quB = " << m_quA << ", " << m_quB << "\n";
    // gsInfo << "MIN JAC G = " << minDJ << "\n";

    if (minDJ <= m_checkForInf_eps) // Check for too small value of detJ at gauss points
    {
        return false;
    }

    return true;
}
