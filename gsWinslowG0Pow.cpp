#include <gismo.h>
#include "gsWinslowG0Pow.h"
using namespace gismo;

real_t gsWinslowG0Pow::evalObj() const {
     //gsInfo << "gsWinslowG0Pow: evalObj() \n";
     //
     //
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(*m_integrationBasis);

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quB",m_quB);
    ev.options().setReal("quA",m_quA);
    A.options().addInt("quRule","quad rule", m_quRule);

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

    geometryMap G0 = A.getMap(m_mp0);

    auto detJ = jac(G).det(); 
    auto detJ0 = jac(G0).det(); 

    real_t q = -2.0/m_mp->domainDim();

    auto J0invJ = jac(G0).inv() * jac(G);

    auto W = (J0invJ)%(J0invJ)*pow(detJ,q)*pow(detJ0,-q+1);

    real_t out = ev.integral(W);

    if (m_addCorners)
    {
        gsMatrix<> corners;
        index_t dim = m_mp->domainDim();
        corners.setZero(dim,pow(2,dim));
        if (dim == 2)
        {
            corners << 0, 0, 1, 1,
                       0, 1, 0, 1;

        } else {
            corners << 0, 0, 0, 0, 1, 1, 1, 1,
                       0, 0, 1, 1, 0, 0, 1, 1,
                       0, 1, 0, 1, 0, 1, 0, 1;
        }
            

        for (index_t i = 0; i < corners.cols(); i++)
        {
            //gsInfo << "i = " << i << " ...  ... ";
            for (index_t p = 0; p < m_mp->nBoxes(); p++)
            {
                //out += m_alpha*ev.eval( jac(G).det(), corners.col(i), p)(0,0);
                out += m_alpha*ev.eval( W , corners.col(i), p)(0,0);
            }
        }
    }


    return out;

}

gsVector<> gsWinslowG0Pow::gradAll(gsDofMapper &space_mapper) const {
    // gsInfo << "gradObj() \n";
    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(*m_mp);
    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);
    A.options().setInt("quB",m_quB);
    A.options().setReal("quA",m_quA);
    A.options().addInt("quRule","quad rule", m_quRule);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(*m_mp);
    geometryMap G0 = A.getMap(m_mp0);

    space u = A.getSpace(dbasis,m_mp->geoDim()); // The gradient is a vector with targetDim values

    A.initSystem();

    real_t q = -2.0/m_mp->domainDim();

    auto detJ = jac(G).det(); 
    auto detJ0 = jac(G0).det(); 

    auto J0invJ = jac(G0).inv() * jac(G);
    auto dJdc = jac(G0).inv() * jac(u);

    auto JTJ = (J0invJ%J0invJ).val();
    auto dJdcTJ = dJdc%J0invJ;

    A.assemble(q*pow(detJ,q)*JTJ*jac(u)%jac(G).inv().tr()*pow(detJ0,-q+1));
    A.assemble(2*pow(detJ,q)*dJdcTJ*pow(detJ0,-q+1));

    auto dWdc = q*pow(detJ,q)*JTJ*jac(u)%jac(G).inv().tr()*pow(detJ0,-q+1) + 2*pow(detJ,q)*dJdcTJ*pow(detJ0,-q+1);

    gsVector<> out = A.rhs();

    if (m_addCorners)
    {
        gsMatrix<> corners;
        index_t dim = m_mp->domainDim();
        corners.setZero(dim,pow(2,dim));
        if (dim == 2)
        {
            corners << 0, 0, 1, 1,
                       0, 1, 0, 1;

        } else {
            corners << 0, 0, 0, 0, 1, 1, 1, 1,
                       0, 0, 1, 1, 0, 0, 1, 1,
                       0, 1, 0, 1, 0, 1, 0, 1;
        }
            

        for (index_t p = 0; p < m_mp->nBoxes(); p++)
        {
                gsMatrix< unsigned > actives;
                dbasis.basis(p).active_into(corners,actives);
                index_t n_actives = actives.rows();

            for (index_t i = 0; i < corners.cols(); i++)
            {
                gsVector<> pt = corners.col(i);

                gsVector<> M = ev.eval( dWdc , pt, p);

                for (index_t d = 0; d < m_mp->targetDim(); d++)
                {
                    gsVector<> Md = M.segment(d*n_actives,n_actives);

                    for(index_t k = 0; k < n_actives; k++)
                    {
                        index_t jj = actives(k,i) + m_patchShift[p] + m_shift_flat[d];

                        out[jj] += m_alpha*Md[k];
                    }

                }

            }
        }
    }

    space_mapper = u.mapper();
    return out;

}
