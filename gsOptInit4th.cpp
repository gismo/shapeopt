#include <gismo.h>
#include <math.h>
#include "gsOptInit4th.h"
using namespace gismo;

gsOptInit4th::gsOptInit4th(gsMultiPatch<>::Ptr mpin, gsMultiPatch<>::Ptr mp_goal):
    gsParamMethod(mpin),
    m_mp_goal(mp_goal)
{
	m_d = m_mp_goal->targetDim();

    n_patches = m_mp_goal->nBoxes();
    n_boundaries = m_mp_goal->nBoundary();
    gsDebugVar(n_boundaries);
    n_design = m_d*(m_d + 1); // N.o entries in A and b in total

    // Start with identity
    m_A.setIdentity(m_d,m_d);
    m_b.setZero(m_d);

    setupOptParameters();
	
};

void gsOptInit4th::setupOptParameters()
{
    // The desing variables is the tagged variables, ie. the boundary controlpoints
    m_numDesignVars = n_design;

	m_curDesign = getDesignVars();

    // Essentially no design bounds
    m_desLowerBounds.setConstant(m_numDesignVars, -1e9);
    m_desUpperBounds.setConstant(m_numDesignVars, 1e9);

	m_numConstraints = 0;
	m_numConJacNonZero = 0;
}

real_t gsOptInit4th::evalObj ( ) const
{
    //gsInfo << "EVALOBJ\n" << std::flush;
    real_t out = 0;

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp_goal->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp_goal->basis(ps.patch).boundary(ps);

		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t i = boundaryDofs[j];

            gsVector<> coef         = m_mp->patch(ps.patch).coef(i);
            gsVector<> coef_goal    = m_mp_goal->patch(ps.patch).coef(i);

			out += 0.5*(m_A*coef + m_b - coef_goal).squaredNorm();
		}
	}

	return out;

}

real_t gsOptInit4th::evalObj ( const gsAsConstVector<real_t> & u) const
{
    //gsInfo << "EVALOBJ\n" << std::flush;
    updateFromDesignVars(u);
    return evalObj();
};

gsVector<> gsOptInit4th::gradObj() const
{
    gsVector<> result(n_design);

    // Set result to zero
    for (index_t kk = 0; kk < n_design; kk++)
        result[kk] = 0;

    // Derivative with respect to A
    for (index_t ii = 0; ii < m_d; ii++)
    {

        for (index_t jj = 0; jj < m_d; jj++)
        {
            gsMatrix<> diffA;
            diffA.setZero(m_d,m_d);

            diffA(ii,jj) = 1;

            index_t kk = getDesignIndex(ii,jj);
             

	        for( index_t k = 0; k < n_boundaries; k++)
	        {
	        	patchSide ps = m_mp_goal->boundaries()[k];
                gsVector<unsigned> boundaryDofs = m_mp_goal->basis(ps.patch).boundary(ps);

	         	for( index_t j = 0; j < boundaryDofs.size(); j++)
	         	{
	         		index_t i = boundaryDofs[j];

                    gsVector<> coef         = m_mp->patch(ps.patch).coef(i);
                    gsVector<> coef_goal    = m_mp_goal->patch(ps.patch).coef(i);

	         		result[kk] += (diffA*coef).transpose()*(m_A*coef + m_b - coef_goal);
                    //result[kk] = 0;
	         	}
	         }
        }
    }

    // Derivative with respect to b
    for (index_t ii = 0; ii < m_d; ii++)
    {

       gsVector<> diffb;
       diffb.setZero(m_d);

       diffb[ii] = 1;

       index_t kk = getDesignIndex(ii);

	   for( index_t k = 0; k < n_boundaries; k++)
	   {
	   	    patchSide ps = m_mp_goal->boundaries()[k];
            gsVector<unsigned> boundaryDofs = m_mp_goal->basis(ps.patch).boundary(ps);

	    	for( index_t j = 0; j < boundaryDofs.size(); j++)
	    	{
	    		index_t i = boundaryDofs[j];

                gsVector<> coef         = m_mp->patch(ps.patch).coef(i);
                gsVector<> coef_goal    = m_mp_goal->patch(ps.patch).coef(i);

	    		result[kk] += diffb.transpose()*(m_A*coef + m_b - coef_goal);
                //result[kk] = 0;
	    	}
	    }
    }
     
    return result;

}


void gsOptInit4th::gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
   // gsInfo << "GRADOBJ\n" << std::flush;
    updateFromDesignVars(u);
    result = gradObj();


};

void gsOptInit4th::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	return;
};

void gsOptInit4th::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	return;
};


void gsOptInit4th::print()
{
  gsInfo << "m_numDesignVars  = " <<  m_numDesignVars << "\n";
  gsInfo << "m_numConstraints  = " <<  m_numConstraints << "\n";

  gsInfo << "m_numConJacNonZero  = " <<  m_numConJacNonZero << "\n";
  gsInfo << "m_conJacRows.size  = " <<  m_conJacRows.size() << "\n";

  // for(std::vector<index_t>::iterator it = m_conJacRows.begin(); it != m_conJacCols.end(); ++it){
  //     gsInfo << " " << *it;
  // }

  gsInfo << "\nm_conJacCols.size  = " <<  m_conJacCols.size() << "\n";

  gsInfo << "m_conUpperBounds.size() = " << m_conUpperBounds.size() << "\n";
  gsInfo << "m_conLowerBounds.size() = " << m_conLowerBounds.size() << "\n";
  gsInfo << "m_desUpperBounds.size() = " << m_desUpperBounds.size() << "\n";
  gsInfo << "m_desLowerBounds.size() = " << m_desLowerBounds.size() << "\n";

  gsInfo << "m_curDesign.size() = " << m_curDesign.size() << "\n";

  gsMatrix<> disp(m_desUpperBounds.size(),2);
  disp << m_desLowerBounds,m_desUpperBounds;
  // gsInfo << ".. design upper and lower bounds\n";
  // gsInfo << disp << "\n";
  //
  gsMatrix<> disp2(m_conUpperBounds.size(),2);
  disp2 << m_conLowerBounds,m_conUpperBounds;
  // gsInfo << ".. constraint upper and lower bounds\n";
  // gsInfo << disp2 << "\n";

}

bool gsOptInit4th::intermediateCallback()
{
    return true;
}

index_t gsOptInit4th::getDesignIndex(const index_t i, const index_t j) const
{
    return i + m_d*j;
};

index_t gsOptInit4th::getDesignIndex(const index_t i) const
{
    return m_d*m_d + i;
};

void gsOptInit4th::updateFromDesignVars( const gsVector<> &u ) const
{
    // Update A
    for (index_t i = 0; i < m_d; i++)
    {
        for (index_t j = 0; j < m_d; j++)
        {
            index_t k = getDesignIndex(i,j);

            m_A(i,j) = u[k];
        }

    }

    // Update b
    for (index_t i = 0; i < m_d; i++)
    {
        index_t k = getDesignIndex(i);

        m_b[i] = u[k];
    }

}

gsVector<> gsOptInit4th::getDesignVars() const
{

    gsVector<> out(n_design);

    // Get A
    for (index_t i = 0; i < m_d; i++)
    {
        for (index_t j = 0; j < m_d; j++)
        {
            index_t k = getDesignIndex(i,j);

            out[k] = m_A(i,j);
        }

    }

    // Update b
    for (index_t i = 0; i < m_d; i++)
    {
        index_t k = getDesignIndex(i);

        out[k] = m_b[i];
    }

    return out;
    
}

void gsOptInit4th::updateMP()
{
    for (index_t p = 0; p < n_patches; p++)
    {
        gsMatrix<> coefs = m_mp->patch(p).coefs();
        for (index_t i = 0; i < coefs.rows(); i++)
        {
            coefs.row(i) = m_A*coefs.row(i).transpose() + m_b;
        }

        m_mp->patch(p).setCoefs(coefs);
    }
            

}
