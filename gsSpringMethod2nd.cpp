#include <gismo.h>
#include "gsSpringMethod2nd.h"
using namespace gismo;

// FIXIT make dimension independent
void gsSpringMethod2nd::setupSystem(gsDofMapper &mapper, gsSparseMatrix<> &A, gsVector<> &b, index_t coord){

	index_t dim = m_mp->targetDim();

    index_t n = mapper.freeSize();

    // Prepare system and rhs
    gsSparseMatrix<> tmp(n,n);
    A.swap(tmp);
    A.reserve( gsVector<int>::Constant(A.cols(), 2*dim) );
    b.setZero(n);

    // For each patch
    for (index_t p = 0; p < m_mp->nBoxes(); p++){
        gsVector<index_t> stride; // Tensor strides
       	gsVector<index_t> size; // Tensor strides

       	gsVector<index_t,2> stride_2; // Tensor strides
       	gsVector<index_t,2> size_2; // Tensor strides
       	gsVector<index_t,3> stride_3; // Tensor strides
       	gsVector<index_t,3> size_3; // Tensor strides

        // Prepare stride and size
		if (dim == 2)
		{
        	static_cast<gsTensorBSplineBasis<2>&> (m_mp->patch(p).basis()).stride_cwise(stride_2);
        	static_cast<gsTensorBSplineBasis<2>&> (m_mp->patch(p).basis()).size_cwise(size_2);
			stride = stride_2;	
			size = size_2;	
		}

		if (dim == 3)
		{
        	static_cast<gsTensorBSplineBasis<3>&> (m_mp->patch(p).basis()).stride_cwise(stride_3);
        	static_cast<gsTensorBSplineBasis<3>&> (m_mp->patch(p).basis()).size_cwise(size_3);
			stride = stride_3;	
			size = size_3;	
		}

        const index_t dd = 2*dim;

        for (int i = 0; i < m_mp->patch(p).coefsSize(); i++)
        {
            index_t ii = mapper.index(i,p);

            if( !mapper.is_free_index(ii))
                continue;

            // if we are on the boundary of the gsMultiPatch we only have (dim-1)*2 neighbors that are 
			// also on the bmnd!!
            // FIXIT: maybe use mappers instead of is_boundary
			std::vector< index_t > dirs;
			getDirections(dirs,i,p);
			index_t n_dirs = dirs.size();

			A(ii,ii) = (dim - n_dirs)*2;

            for ( unsigned k = 0; k<dim; k++ ) // for all neighbors
            {
                for ( int s = -1; s<2; s+=2 ) // +/- (up or down)
                {
					bool stop = false;
					for (index_t t = 0; t < n_dirs; t++)
					{
						if (k == dirs[t]) stop = true;

					}
					if (stop) continue;

                    const unsigned j = i + s * stride[k];

                    //if (j < 0 || j >= m_mp->patch(p).coefsSize()) // Out of bounds?
                    //{
                    //    continue;
                    //}

                    index_t jj = mapper.index(j,p);// get j-index in the matrix

                    if ( mapper.is_free_index(jj) )
                    {
                       	A(ii,jj) = -1;
						
                        // gsInfo << ii << ", " << jj << " in A is altered...\n";
                    }
                    else // boundary node
                    {
                        // for 2D: add only half on the interfaces where this is traversed again later..
                        // b[ii] += (mapper.is_coupled_index(ii) ? 0.5 : 1.0) * m_mp->patch(p).coef(j,coord);

                        // FIXIT think about whether this strategy is correct or not
                        if (mapper.is_coupled_index(ii)) { // If we are on interface
                            // gsInfo << ii << " is coupled to ";
                       
                            // Count no. of times the neighbor is coupled, this is the number of time we will visit it
                            std::vector< std::pair<index_t, index_t > > vec;
                            mapper.preImage(jj,vec);
                            real_t n_coupled = vec.size();
                       
                            // Only add 1/n_coupled of the coef since it will be visited n_coupled times.
                            b[ii] += 1.0/n_coupled * m_mp->patch(p).coef(j,coord);
                       
                            // gsInfo << vec.size() << "\n";
                        } else { // If we are not on an interface, we can just add the neighbor to the rhs.
                            b[ii] +=  m_mp->patch(p).coef(j,coord);
                        }
                    }
                }
            }
        }
    }
	std::string name = "A";
    std::ofstream file (name);
    for(index_t i = 0; i < A.rows(); i++){
       for (index_t j = 0; j < A.cols(); j++){
           file << A(i,j);
           file << " ";
       }
       file << std::setprecision(12) << "\n";
    }
    file.close();

}

void gsSpringMethod2nd::getDirections(std::vector<index_t> &result, index_t i, index_t p)
{

    for (index_t k = 0; k < m_mp->nBoundary(); k++){
        patchSide ps = m_mp->boundaries()[k];
        if (ps.patch != p) continue;

        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
        for (index_t j = 0; j < boundaryDofs.size(); j++){
            if (boundaryDofs[j] == i) result.push_back(ps.direction());
        }
    }

}

bool gsSpringMethod2nd::is_double_boundary(index_t i, index_t p)
{
	index_t count = 0;
    for (index_t k = 0; k < m_mp->nBoundary(); k++){
        patchSide ps = m_mp->boundaries()[k];
        if (ps.patch != p) continue;

        gsVector<unsigned> boundaryDofs = m_mp->basis(ps.patch).boundary(ps);
        for (index_t j = 0; j < boundaryDofs.size(); j++){
            if (boundaryDofs[j] == i) count++;
        }
    }
    return count > 1;
}
