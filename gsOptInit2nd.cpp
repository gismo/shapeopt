#include <gismo.h>
#include <math.h>
#include "gsOptInit2nd.h"
using namespace gismo;

gsOptInit2nd::gsOptInit2nd(memory::shared_ptr<gsMultiPatch<>> mpin):
    gsParamMethod(mpin)
{
	n_boundaries = m_mp->nBoundary();
	m_d = m_mp->targetDim();
	n_corners = pow(2,m_d);

    setupOptParameters();
	
	getCorners();

};

void gsOptInit2nd::setupOptParameters()
{
    // The desing variables is the tagged variables, ie. the boundary controlpoints
    m_numDesignVars = 2*m_d;
	gsVector<> zers;
	zers.setZero(m_numDesignVars);

    m_curDesign = zers;

    // Essentially no design bounds
    m_desLowerBounds.setConstant(m_numDesignVars, -1e9);
    m_desUpperBounds.setConstant(m_numDesignVars, 1e9);

	m_numConstraints = 0;
}

real_t gsOptInit2nd::evalObj ( const gsAsConstVector<real_t> & u) const
{
	real_t out = 0;

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		index_t l = m_direction_map[ps.direction()];
		
		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t i = boundaryDofs[j];
			index_t ii = m_mappers[l].index(i,ps.patch);
			if (isCorner(ii))
			{
				//gsDebugVar(i);
				//gsDebugVar(k);
				//gsDebugVar(l);
				//gsDebugVar(m_mp->patch(ps.patch).coef(i,l));
				out += 0.5*pow(u[k] - m_mp->patch(ps.patch).coef(i,l),2);
			}
		}
	}

	return out;
};


void gsOptInit2nd::gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	gsVector<> out;
	out.setZero(m_numDesignVars);

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		index_t l = m_direction_map[ps.direction()];
		
		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t i = boundaryDofs[j];
			index_t ii = m_mappers[l].index(i,ps.patch);
			if (isCorner(ii))
				out[k] += u[k] - m_mp->patch(ps.patch).coef(i,l);
		}
	}

	result = out;

};

void gsOptInit2nd::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	return;
};

void gsOptInit2nd::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	return;
};


void gsOptInit2nd::print()
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

bool gsOptInit2nd::intermediateCallback()
{
    return true;
}

void gsOptInit2nd::getCorners()
{
	// Count to find which cps is corners
	m_count.setZero(n_cps/m_d);

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t ii = m_mappers[0].index(boundaryDofs[j],ps.patch);
			m_count[ii]++;
		}
	}

	// Create direction map for each side.
	m_direction_map.setZero(m_d);
	std::vector< bool > para_direction_flag(m_d);
	std::vector< bool > phys_direction_flag(m_d);

 	for (index_t d = 0; d < m_d; d++)
	{
		para_direction_flag[d] = false;
		phys_direction_flag[d] = false;
	}

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		index_t l = ps.direction();
		if (para_direction_flag[l]) // Stop if this direction is already considered.
			continue;

		gsVector<> avg;
		avg.setZero(m_d);
		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t i = boundaryDofs[j];
			index_t ii = m_mappers[0].index(i,ps.patch);
			if(isCorner(ii))
			{
				for(index_t d = 0; d < m_d; d++)
				{
					avg[d] += m_mp->patch(ps.patch).coef(i,d)/pow(2,m_d-1);
				}
			}
		}

		gsVector<> var;
		var.setZero(m_d);

		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t i = boundaryDofs[j];
			index_t ii = m_mappers[0].index(i,ps.patch);
			if(isCorner(ii))
			{
				// Calculate variances
				for(index_t d = 0; d < m_d; d++)
				{
					var[d] += pow(m_mp->patch(ps.patch).coef(i,d) - avg[d],2);
				}
				
				
			}
		}

		// Find direction with minimal variance
		real_t min = 1e6;
	
		gsDebugVar(l);
		for(index_t d = 0; d < m_d; d++)
		{
			gsDebugVar(d);
			gsDebugVar(var[d]);
			gsDebugVar(avg[d]);
			if (phys_direction_flag[d]) continue;

			if (var[d] < min)
			{
				min = var[d];
				m_direction_map[l] = d;
			}
		}
		phys_direction_flag[m_direction_map[l]] = true;
		para_direction_flag[l] = true;
	}

 	for (index_t d = 0; d < m_d; d++)
	{
		gsInfo << "m_direction_map[ " << d << " ] = " << m_direction_map[d] << "\n"; 
	}

	// Map each boundary to a design variable
	m_boundary_to_dof.setZero(n_boundaries);

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		index_t l = m_direction_map[ps.direction()];
	}
}

void gsOptInit2nd::updateCorners(gsMultiPatch<>::Ptr mp)
{
	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		index_t l = m_direction_map[ps.direction()];
		
		for( index_t j = 0; j < boundaryDofs.size(); j++)
		{
			index_t i = boundaryDofs[j];
			index_t ii = m_mappers[0].index(i,ps.patch);
			if (isCorner(ii))
				mp->patch(ps.patch).coef(i,l) = m_curDesign(k,0);
		}
	}
}

bool gsOptInit2nd::isCorner(index_t ii) const
{
	return m_count[ii] >= m_d && !m_mappers[0].is_coupled(ii);
}
