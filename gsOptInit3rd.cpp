#include <gismo.h>
#include <math.h>
#include "gsOptInit3rd.h"
using namespace gismo;

gsOptInit3rd::gsOptInit3rd(memory::shared_ptr<gsMultiPatch<>> mpin):
    gsParamMethod(mpin)
{
	n_boundaries = m_mp->nBoundary();
	m_d = m_mp->targetDim();
	n_corners = pow(2,m_d);
	n_boxes = m_mp->nBoxes();

    setupOptParameters();
	
	getCorners();

};

void gsOptInit3rd::setupOptParameters()
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

real_t gsOptInit3rd::evalObj ( const gsAsConstVector<real_t> & u) const
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
				index_t dof = m_boundary_to_dof(ps.patch,ps.index());
				out += 0.5*pow(u[dof] - m_mp->patch(ps.patch).coef(i,l),2);
			}
		}
	}

	return out;
};


void gsOptInit3rd::gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
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
			{
				index_t dof = m_boundary_to_dof(ps.patch,ps.index());
				out[dof] += u[dof] - m_mp->patch(ps.patch).coef(i,l);
			}
		}
	}

	result = out;

};

void gsOptInit3rd::evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	return;
};

void gsOptInit3rd::jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
{
	return;
};


void gsOptInit3rd::print()
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

bool gsOptInit3rd::intermediateCallback()
{
    return true;
}

void gsOptInit3rd::getCorners()
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
	
		for(index_t d = 0; d < m_d; d++)
		{
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

	// Map each boundary to a design variable
	m_boundary_to_dof.setZero(n_boxes,2*m_d+1);
	gsMatrix< bool > ps_flag(n_boxes,2*m_d+1);

	for( index_t p = 0; p < n_boxes; p++)
		for ( index_t i = 0; i < 2*m_d + 1; i++)
			ps_flag(p,i) = false;

	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);

		if (ps_flag(ps.patch,ps.index()))
			continue;

		markBoundaryAndOppositeToDof(ps,ps_flag);

		std::vector< patchSide > neighbours = getUnmarkedNeighbours(ps,ps_flag);

		for( std::vector< patchSide >::iterator it = neighbours.begin(); it != neighbours.end(); ++it)
			markBoundaryAndOppositeToDof(*it,ps_flag);
	}

 	for (index_t d = 0; d < m_d; d++)
	{
		gsInfo << "m_direction_map[ " << d << " ] = " << m_direction_map[d] << "\n"; 
	}

	
	for( index_t k = 0; k < n_boundaries; k++)
	{
		patchSide ps = m_mp->boundaries()[k];
		gsInfo << "patch " << ps.patch << " index " << ps.index() << " is mapped to " << m_boundary_to_dof(ps.patch,ps.index()) << "\n";

	}

}

void gsOptInit3rd::updateCorners(gsMultiPatch<>::Ptr mp)
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
			{	
				index_t dof = m_boundary_to_dof(ps.patch,ps.index());
				mp->patch(ps.patch).coef(i,l) = m_curDesign(dof,0);
			}
		}
	}
}

bool gsOptInit3rd::isCorner(index_t ii) const
{
	return m_count[ii] >= m_d && !m_mappers[0].is_coupled(ii);
}

void gsOptInit3rd::markBoundaryAndOppositeToDof(patchSide &ps, gsMatrix< bool > &ps_flag)
{
	index_t p = ps.patch;
	index_t l = m_direction_map[ps.direction()];

	// Mark this boundary to k
	markBoundaryToDof(ps, p, l, ps_flag);

	boxSide opposite = ps.opposite();

	markBoundaryToDof(opposite, p, l + m_d, ps_flag);
}

void gsOptInit3rd::markBoundaryToDof(boxSide &bs, index_t p , index_t dof, gsMatrix< bool > &ps_flag)
{
	index_t i = bs.index();
	m_boundary_to_dof(p, i) = dof;
	ps_flag(p, i) = true;
}

std::vector< patchSide > gsOptInit3rd::getUnmarkedNeighbours(patchSide &ps, gsMatrix< bool > &ps_flag)
{
	std::vector< patchSide > out;

    gsMatrix<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
	gsMatrix<unsigned> globalDofs(boundaryDofs.size(),1);

	m_mappers[0].localToGlobal(boundaryDofs,ps.patch,globalDofs);

	sort(globalDofs);

	for( index_t k = 0; k < n_boundaries; k++ )
	{
		patchSide ps = m_mp->boundaries()[k];

		if (ps_flag(ps.patch,ps.index()) || ps.patch == k)
			continue;

        gsMatrix<unsigned> boundaryDofs_k = m_mp->basis(ps.patch).boundary(ps);
		gsMatrix<unsigned> globalDofs_k(boundaryDofs.size(),1);

		m_mappers[0].localToGlobal(boundaryDofs_k,ps.patch,globalDofs_k);

		sort(globalDofs_k);

		index_t count = countIntersection(globalDofs_k,globalDofs);

		if (count >= m_d-1)
		{
			out.push_back(ps);
			ps_flag(ps.patch,ps.index()) = true;
		}
	}

	return out;
}

void gsOptInit3rd::sort(gsMatrix< unsigned > &v)
{
	std::sort(v.data(),v.data()+v.size());
}

index_t gsOptInit3rd::countIntersection(gsMatrix< unsigned > &v, gsMatrix< unsigned > &w)
{
	index_t count = 0;
	index_t k = 0;
	for (index_t i = 0; i < v.size(); i++)
	{
		gsDebugVar(i);
		gsDebugVar(k);
		while(k < w.size() && v(i,0) < w(k,0))
		{
			k++;
			if (k >= w.size())
				return count;
		}

		if (v(i,0) == w(k,0))
			count++;
		
	}
	
	return count;
}

