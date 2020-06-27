#include <gismo.h>
#include "gsWinslowPow.h"
using namespace gismo;

real_t gsWinslowPow::evalObj() const {
     //gsInfo << "gsWinslowPow: evalObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(*m_integrationBasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);

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

    auto detJ = jac(G).det(); 

    real_t q = -2.0/m_mp->domainDim();

    real_t out = ev.integral(jac(G)%jac(G)*pow(detJ,q));

    //gsDebugVar(out);

    //real_t outNoPow = ev.integral(jac(G)%jac(G)*pow(detJ,-1));
    //gsDebugVar(outNoPow);

    return out;

}

gsVector<> gsWinslowPow::gradAll(gsDofMapper &space_mapper) const {
    // gsInfo << "gradObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(*m_integrationBasis);

    gsExprEvaluator<> ev(A);
    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    space u = A.getSpace(dbasis,m_mp->geoDim()); // The gradient is a vector with targetDim values

    A.initSystem();

    real_t q = -2.0/m_mp->domainDim();

    auto detJ = jac(G).det(); 
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
    A.assemble(q*pow(detJ,q)*JTJ*jac(u)%jac(G).inv().tr());
    A.assemble(2*pow(detJ,q)*jac(u)%jac(G));


    space_mapper = u.mapper();
    return A.rhs();

}

gsMatrix<> gsWinslowPow::hessAll(gsDofMapper &space_mapper) const {

        gsExprAssembler<> A(1,1);
        gsMultiBasis<> dbasis(*m_mp);
        A.setIntegrationElements(*m_integrationBasis);
        gsExprEvaluator<> ev(A);

        A.options().setInt("quB",m_quB);
        A.options().setReal("quA",m_quA);

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;

        geometryMap G = A.getMap(*m_mp);

        space u = A.getSpace(dbasis,m_mp->geoDim());

        A.initSystem();

        real_t q = -2.0/m_mp->domainDim();

        auto detJ= jac(G).det(); // The inverse of abs(det J)JTJ

        auto JTJ = (jac(G)%jac(G)).val();

        auto W = JTJ*pow(detJ,q);
        auto dWdc = 2*pow(detJ,q)*jac(u)%jac(G) 
                      + q*pow(detJ,q)*JTJ*jac(u)%jac(G).inv().tr();

        auto dJinvTdc = - matrix_by_space_tr(jac(G).inv(),jac(u))*jac(G).inv().tr();

        A.assemble(2*pow(detJ,q)*jac(u)%jac(u));
        //
        A.assemble(2*q*pow(detJ,q)*(jac(u)%jac(G))*(jac(u)%jac(G).inv().tr()).tr());
        A.assemble(2*q*pow(detJ,q)*(jac(u)%jac(G).inv().tr())*(jac(u)%jac(G)).tr());
        // A.assemble(-2*(jac(u)%jac(G).inv().tr())*(jac(G)%jac(u)).tr()*detJinv);
        //
        A.assemble(q*q*W*(jac(u)%jac(G).inv().tr())*(jac(u)%jac(G).inv().tr()).tr());
        //
        // // A.assemble(-detJinv2*JTJ*jac(u)%jac(u).adj().tr());
        A.assemble(-q*W*jac(u).tr()%(jac(G).inv()*jac(u)*jac(G).inv()));

        space_mapper = u.mapper();
        return A.matrix();
    }


