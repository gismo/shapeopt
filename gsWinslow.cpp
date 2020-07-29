#include <gismo.h>
#include "gsWinslow.h"
using namespace gismo;

gsWinslow::gsWinslow(memory::shared_ptr<gsMultiPatch<>> mpin, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps):
    gsOptParamMethod(mpin,use_dJC,useTensorStructureforDJC),
    m_checkForInf(checkForInf),
    m_checkForInf_eps(checkForInf_eps)
{
    //gsInfo << "Check for inf in Winslow: " << checkForInf << " with eps " << checkForInf_eps << "\n";
}

gsWinslow::gsWinslow(memory::shared_ptr<gsMultiPatch<>> mpin, std::vector< gsDofMapper > mappers, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps):
    gsOptParamMethod(mpin,mappers,use_dJC,useTensorStructureforDJC),
    m_checkForInf(checkForInf),
    m_checkForInf_eps(checkForInf_eps)
{
    //gsInfo << "Check for inf in Winslow: " << checkForInf << " with eps " << checkForInf_eps << "\n";
}

real_t gsWinslow::evalObj() const {
    // gsInfo << "evalObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(*m_integrationBasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);

    ev.options().addInt("quRule","quad rule", m_quRule);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    real_t minDJ = ev.min(jac(G).det());

    if (m_checkForInf)
    {
        if (minDJ <= m_checkForInf_eps) // Check for too small value of detJ at gauss points
        {
            // gsWriteParaview(*m_mp,"geomwin");
            // gsInfo << "norm of flat : " << getFlat().norm() << "\n";
            // gsInfo << "norm of tagged : " << getTagged().norm() << "\n";
            // gsInfo << "quA, quB : " << m_quA << ", " << m_quB << "\n";
            // gsInfo << "Min detJ in winslow: " << minDJ << "\n\n";
            return std::numeric_limits<double>::infinity();
        }
    }

    auto detJinv = jac(G).inv().det(); // The inverse of det J

    real_t out = ev.integral(jac(G)%jac(G)*detJinv);

    return out;

}

real_t gsWinslow::evalObj(index_t p) const {
    // gsInfo << "evalObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(m_mp->patch(p));
    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);
    ev.options().addInt("quRule","quad rule", m_quRule);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    real_t minDJ = ev.min(jac(G).det());

    if (m_checkForInf)
    {
        if (minDJ <= m_checkForInf_eps) // Check for too small value of detJ at gauss points
        {
            // gsWriteParaview(*m_mp,"geomwin");
            // gsInfo << "norm of flat : " << getFlat().norm() << "\n";
            // gsInfo << "norm of tagged : " << getTagged().norm() << "\n";
            // gsInfo << "quA, quB : " << m_quA << ", " << m_quB << "\n";
            // gsInfo << "Min detJ in winslow: " << minDJ << "\n\n";
            return std::numeric_limits<double>::infinity();
        }
    }

    auto detJinv = jac(G).inv().det(); // The inverse of det J

    real_t out = ev.integral(jac(G)%jac(G)*detJinv);

    return out;

}

gsVector<> gsWinslow::gradObj() const{
    gsDofMapper mapper;
    gsVector<> all = gradAll(mapper);
    return mapMatrix(mapper,all);

}

gsMatrix<> gsWinslow::hessAll(gsDofMapper &space_mapper) const
{
        gsExprAssembler<> A(1,1);
        gsMultiBasis<> dbasis(*m_mp);
        A.setIntegrationElements(*m_integrationBasis);
        gsExprEvaluator<> ev(A);

        A.options().setInt("quB",m_quB);
        A.options().setReal("quA",m_quA);
        A.options().addInt("quRule","quad rule", m_quRule);

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;

        geometryMap G = A.getMap(*m_mp);

        space u = A.getSpace(dbasis,m_mp->geoDim());

        A.initSystem();

        auto detJinv = jac(G).inv().det().val(); // The inverse of abs(det J)JTJ

        auto JTJ = (jac(G)%jac(G)).val();

        auto W = JTJ*detJinv;
        auto dWdc = 2*detJinv*(jac(u)%jac(G))
                        - detJinv*(jac(u)%jac(G).inv().tr())*JTJ;

        auto dJinvTdc = - matrix_by_space_tr(jac(G).inv(),jac(u))*jac(G).inv().tr();

        // A.assemble(2*detJinv*jac(u)%jac(u));
        // A.assemble(-2*(jac(u)%jac(G).inv().tr())*(jac(u)%jac(G)).tr()*detJinv);

        // A.assemble(- jac(u) % (W.val() * dJinvTdc));

        // A.assemble(- (jac(u)%jac(G).inv().tr()) * dWdc.tr());
        // A.assemble(- 2*detJinv*(jac(u)%jac(G).inv().tr())*(jac(u)%jac(G)).tr());
        // A.assemble(W*(jac(u)%jac(G).inv().tr()) * (jac(u)%jac(G).inv().tr()).tr());

        //========================//

        A.assemble(2*detJinv*jac(u)%jac(u));
        //
        A.assemble(-2*detJinv*(jac(u)%jac(G))*(jac(u)%jac(G).inv().tr()).tr());
        A.assemble(-2*detJinv*(jac(u)%jac(G).inv().tr())*(jac(u)%jac(G)).tr());
        // A.assemble(-2*(jac(u)%jac(G).inv().tr())*(jac(G)%jac(u)).tr()*detJinv);
        //
        A.assemble(W*(jac(u)%jac(G).inv().tr())*(jac(u)%jac(G).inv().tr()).tr());
        //
        // // A.assemble(-detJinv2*JTJ*jac(u)%jac(u).adj().tr());
        A.assemble(W*jac(u).tr()%(jac(G).inv()*jac(u)*jac(G).inv()));

        space_mapper = u.mapper();
        return A.matrix();
    }

gsMatrix<> gsWinslow::hessAll() const{
        gsDofMapper tmp;
        return hessAll(tmp);
    }

gsVector<> gsWinslow::gradAll() const{
        gsDofMapper tmp;
        return gradAll(tmp);
    }

gsVector<> gsWinslow::gradAll(gsDofMapper &space_mapper) const{
    // gsInfo << "gradObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(*m_integrationBasis);

    gsExprEvaluator<> ev(A);
    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);
    A.options().addInt("quRule","quad rule", m_quRule);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    space u = A.getSpace(dbasis,m_mp->geoDim()); // The gradient is a vector with targetDim values

    A.initSystem();

    auto detJinv = jac(G).inv().det(); // The inverse of abs(det J)JTJ
    auto detJinv2 = detJinv*detJinv; // The inverse of abs(det J)JTJ
    auto JTJ = (jac(G)%jac(G)).val();

    // A.assemble(2*matrix_by_space(jac(G).tr(),jac(u)).trace()*detJinv
    //     - matrix_by_space(jac(G).inv(),jac(u)).trace()*
    //     (jac(G).tr()*jac(G)).trace().val()*detJinv);
    // A.assemble
    // (
    //      2*detJinv*(jac(u)%jac(G))
    //      - detJinv*(jac(u)%jac(G).inv().tr())*JTJ
    // );
    // A.assemble(2*detJinv*jac(u)%jac(G) - detJinv2*JTJ*jac(u)%jac(G).adj().tr());
    A.assemble(2*detJinv*jac(u)%jac(G) - detJinv*JTJ*jac(u)%jac(G).inv().tr());


    space_mapper = u.mapper();
    return A.rhs();

}

gsVector<> gsWinslow::gradObj(gsVector<> &gradObjTagged) const{

    gsDofMapper mapper;
    gsVector<> all = gradAll(mapper);

    // Map tagged part
    gradObjTagged = mapMatrixToTagged(mapper,all);

    return mapMatrix(mapper,all);

}

real_t gsWinslow::minDetJInGaussPts(index_t incPts){

    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB + incPts);
    ev.options().setReal("quA",m_quA);
    ev.options().addInt("quRule","quad rule", m_quRule);

    typedef gsExprAssembler<>::geometryMap geometryMap;

    geometryMap G = A.getMap(*m_mp);
	//gsDebugVar(ev.max((jac(G).det())));
	//gsDebugVar(ev.min((jac(G).det())));

    return ev.min(jac(G).det());
}

real_t gsWinslow::maxDetJInGaussPts(index_t incPts){

    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB + incPts);
    ev.options().setReal("quA",m_quA);
    ev.options().addInt("quRule","quad rule", m_quRule);

    typedef gsExprAssembler<>::geometryMap geometryMap;

    geometryMap G = A.getMap(*m_mp);
	//gsDebugVar(ev.max((jac(G).det())));
	//gsDebugVar(ev.min((jac(G).det())));

    return ev.max(jac(G).det());
}

void gsWinslow::computeWinslowPerPatch()
{
    m_winslow_per_patch.setZero(m_mp->nBoxes());

    for (index_t p = 0; p < m_mp->nBoxes(); p++)
    {
        m_winslow_per_patch[p] = evalObj(p);
    }
}

