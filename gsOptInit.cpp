#include <gismo.h>
#include "gsOptInit.h"
using namespace gismo;

gsOptInit::gsOptInit(memory::shared_ptr<gsMultiPatch<>> mpin, gsVector<> &tagged_goal):
    gsParamMethod(mpin), m_tagged_goal(tagged_goal)
{
	n_boundaries = m_mp->nBoundary();
	n_tagged_cps = m_mappers[0].taggedSize(); 

	m_global_to_tagged.setConstant(n_cps,-1);
	index_t ind = 0;
    const std::vector< index_t > t = m_mappers[0].getTagged(); // Get tagged from mapper
    for (std::vector< index_t >::const_iterator it = t.begin(); it != t.end(); ++it)
    {
        // *it returns the global index from the mapper
        m_global_to_tagged[*it] = ind;
        ind++;
    }

    setupOptParameters();
	

};

void gsOptInit::setupOptParameters()
{
    // The desing variables is the tagged variables, ie. the boundary controlpoints
    m_numDesignVars = n_tagged + n_boundaries;
	gsVector<> zers;
	zers.setZero(m_numDesignVars);
	zers.segment(0,n_tagged) = getTagged();

    m_curDesign = zers;

    // Essentially no design bounds
    m_desLowerBounds.setConstant(m_numDesignVars, -1e9);
    m_desUpperBounds.setConstant(m_numDesignVars, 1e9);

	for(index_t i = 0; i < n_boundaries; i++)
	{
		//m_desLowerBounds[n_tagged + i] = 2;
		//m_curDesign(n_tagged + i,0) = 2;
	}

	m_numConstraints = 0;
	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		m_numConstraints += boundaryDofs.size();
	}

	// ...that equals 0
    m_conLowerBounds.setZero(m_numConstraints);
    m_conUpperBounds.setZero(m_numConstraints);

	computeConstraintMatrix();

    // compute jac structure
    computeJacStructure();
}

real_t gsOptInit::evalObj ( const gsAsConstVector<real_t> & u) const
{
    // gsInfo << "evalObj\n" << std::flush;
	updateTagged(u);
    return 0.5*(getTagged() - m_tagged_goal).squaredNorm();
};

gsIpOptSparseMatrix gsOptInit::jacobCon() const
{
	return *m_A_ipopt;
}

void gsOptInit::gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "gradObj_into\n" << std::flush;
	updateTagged(u);
	gsDebugVar(u.segment(n_tagged,n_boundaries));

	gsVector<> out;
	out.setZero(m_numDesignVars);

    out.segment(0,n_tagged) = getTagged() - m_tagged_goal;

	result = out;
};

void gsOptInit::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	result = m_A*u;
};

void gsOptInit::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
    // gsInfo << "jacobCon_into\n" << std::flush;
    updateTagged(u);
    gsIpOptSparseMatrix J = jacobCon();
    result = J.values();
};

void gsOptInit::computeJacStructure()
{
    // Use the sparsity provided by gsDetJacConstraint m_dJC
    gsIpOptSparseMatrix J = jacobCon();

    m_numConJacNonZero = J.nnz();
    m_conJacRows = J.rows();
    m_conJacCols = J.cols();

};

void gsOptInit::print()
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

bool gsOptInit::intermediateCallback()
{
    return true;
}

void gsOptInit::computeConstraintMatrix()
{
	m_A.setZero(m_numConstraints,m_numDesignVars);

	index_t slack_ind = 0;
	index_t const_ind = 0;

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		index_t l = ps.direction();
		
		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t jj = m_mappers[l].index(boundaryDofs[j],ps.patch);

			GISMO_ASSERT(m_global_to_tagged[jj] != -1,"WOOOOPS");

			index_t ii = m_global_to_tagged[jj] + m_shift_tagged[l];
			gsInfo << "m_shifttagged " << m_shift_tagged[l] << " is fixed in direction " << l << "\n";
			gsDebugVar(ii);

			gsVector<> tagged = getTagged();
			tagged[ii] = -5;
			updateTagged(tagged);
			m_curDesign.block(0,0,n_tagged,1) = tagged;
			
			m_A(const_ind,ii) = 1;
			m_A(const_ind,slack_ind + n_tagged) = -1;
			const_ind++;
		}
		slack_ind++;
		break;
	}

	m_A_ipopt = memory::make_unique(new gsIpOptSparseMatrix(m_A,1e-1)); 

}
