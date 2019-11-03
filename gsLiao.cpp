#include <gismo.h>
#include "gsLiao.h"
using namespace gismo;

real_t gsLiao::evalObj() const {
    // gsInfo << "evalObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);
    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);

    return ev.integral(fform(G).sqNorm());
}

gsVector<> gsLiao::gradObj() const{
    // gsInfo << "gradObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

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

    A.assemble
    (
         4*((jac(G).tr()*jac(u)) % fform(G))
    );

    gsVector<> all = A.rhs();
    return mapMatrix(u.mapper(),all);

}

gsMatrix<> gsLiao::hessAll(gsDofMapper &space_mapper) const
{
        gsExprAssembler<> A(1,1);
        gsMultiBasis<> dbasis(*m_mp);
        A.setIntegrationElements(dbasis);
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

        // 4*((jac(G).tr()*jac(u)) % fform(G))

        A.assemble
        (
             4*
             (
                 (jac(u)*jac(G).tr()) % (jac(u)*jac(G).tr())
                 +
                 (jac(G).tr()*jac(u)) % (jac(G).tr()*jac(u))
                 +
                 (jac(G).tr()*jac(u)) % (jac(u).tr()*jac(G))
             )
             // jac(u) % jac(u).tr()
        );

        space_mapper = u.mapper();
        return A.matrix();
    }
