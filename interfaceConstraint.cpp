#include <gismo.h>
#include "interfaceConstraint.h"
#include "IpOptSparseMatrix.h"
using namespace gismo;

interfaceConstraint::interfaceConstraint(gsMultiPatch<>* mpin){
    mp = mpin;
    aBigNumber = 100000;

    // Count the total number of controlpoints
      n_controlpoints = 0;
      for(int i = 0; i < mp->nBoxes(); i++){
        n_controlpoints += mp->patch(i).coefsSize();
        // gsInfo << "Patch " << i << " has " << mp->patch(i).coefsSize() << " cps... \n";
      }
      // gsInfo << "n_controlpoints = " << n_controlpoints << "\n \n";

    gsMatrix<> mat = generateConstraintMatrix();
    m_constraintMatrix.swap(IpOptSparseMatrix(mat));
    bndFixInfo.setZero(mp->nBoundary()); // Initialize bndFixInfo to 0, i.e. all boundaries are fixed by default
  };

gsMatrix<> interfaceConstraint::generateConstraintMatrix(){
  // Count the number of constraints
  n_constraints = 0;
  for(int i = 0; i < mp->nInterfaces(); i++){
    boundaryInterface interface = mp->bInterface(i);
    patchSide ps = interface.first();
    n_constraints += 2*(*mp->patch(ps.patch).boundary(ps)).coefsSize();
  }
  // gsInfo << "Number of constraints : " << n_constraints << "\n";
  gsMatrix<real_t> interfaceConstraintMatrix(n_constraints,n_controlpoints*2);
  interfaceConstraintMatrix.setZero(n_constraints,n_controlpoints*2);

  index_t count = 0;
  for(int i = 0; i < mp->nInterfaces(); i++){
    boundaryInterface interface = mp->bInterface(i);
    gsVector<index_t> I = getPatchSideIndicies(interface.first());
    gsVector<index_t> J = getPatchSideIndicies(interface.second());

    // Change direction corresponding to the directionMap
    int firstDirection = 0;
    //If interface at patch 0 is east or west indicate that firstDirection is the second direction
    if (interface.first().index() < 3){
      firstDirection++;
    }
    if (not interface.dirOrientation()[firstDirection]){
      gsVector<index_t> tmp(J);
      for(index_t k = 0; k < J.size(); k++){
          J[k] = tmp[J.size()-k-1];   // Reverse the order
      }
    }

    // gsInfo << "i = " << i << "\n";
    for(index_t k = 0; k < I.size(); k++){
      index_t i = I[k];
      index_t j = J[k];
      addInterfaceConstraint(i,j,count,interfaceConstraintMatrix);
      addInterfaceConstraint(i+n_controlpoints,j+n_controlpoints,count,interfaceConstraintMatrix);

    }
  }
  return interfaceConstraintMatrix;
}

void interfaceConstraint::addInterfaceConstraint(index_t i, index_t j, index_t &count, gsMatrix<> &interfaceConstraintMatrix){
      // Add constraint line at row "count" when to points are connected
      // gsInfo << " ..... " << "Constraint no. " << count <<" : " << i << " <-> " << j << "\n";
      interfaceConstraintMatrix(count,i) = 1;
      interfaceConstraintMatrix(count,j) = -1;
      count++;
}

gsVector<index_t> interfaceConstraint::getPatchSideIndicies(patchSide ps){
  // Generate offset for the given patch
  int n_out = (*mp->patch(ps.patch).boundary(ps)).coefsSize();
  size_t n_coefs = mp->patch(ps.patch).coefsSize();
  gsVector<index_t> patchIndicies(n_out);
  index_t patchOffset = 0;
  for(int i = 0; i < ps.patch; i++){
    patchOffset += mp->patch(i).coefsSize();
  }

  getLocalPatchSideIndicies(ps,patchIndicies);

  for(int i = 0; i < n_out; i++){
    patchIndicies[i] += patchOffset;
  }
  return patchIndicies;
}

void interfaceConstraint::getLocalPatchSideIndicies(patchSide ps, gsVector<index_t> &patchIndicies){
  size_t n = (*mp->patch(ps.patch).boundary(ps)).coefsSize();   // Number of basis functions at edge
  size_t m = mp->patch(ps.patch).coefsSize()/n;                 // Number of basis function at other direction
  int index = ps.index();

  int j = 0;

  switch (index){
    case 1: for(index_t i = 0; i < n*m; i += m){
              patchIndicies[j] = i;
              j++;
            }
            break;
    case 2: for(index_t i = m-1; i < n*m; i += m){
            patchIndicies[j] = i;
            j++;
            }
            break;
    case 3: for(index_t i = 0; i < n; i++){
              patchIndicies[j] = i;
              j++;
            }
            break;
    case 4: for(index_t i = (m-1)*n; i < n*m; i++){
              patchIndicies[j] = i;
              j++;
            }
            break;
  }
}

gsVector<> interfaceConstraint::getUpperBounds(){
  gsVector<> upperBounds;
  upperBounds.setOnes(2*n_controlpoints);
  upperBounds *= aBigNumber;

  // Fix all boundaries of the multipatch domain
  for(index_t i = 0; i < mp->nBoundary(); i++){
    patchSide ps = mp->boundaries()[i];
    getBounds(ps,upperBounds,bndFixInfo[i]);
  }
  return upperBounds;
}

gsVector<> interfaceConstraint::getUpperBounds(index_t fixedPatch){

  gsVector<> upperBounds = getUpperBounds();

  gsMultiPatch<> fixedGeom(mp->patch(fixedPatch));

  for(index_t i = 0; i < fixedGeom.nBoundary(); i++){
    patchSide ps = fixedGeom.boundaries()[i];
    ps.patch = fixedPatch;
    getBounds(ps,upperBounds,0); // 0 means fix in both y and y direction
  }

  return upperBounds;
}

gsVector<> interfaceConstraint::getLowerBounds(gsVector<> upperBounds){
  gsVector<> lowerBounds(upperBounds);
  for(index_t i = 0; i < lowerBounds.size();i++){
    if (lowerBounds[i] == aBigNumber) lowerBounds[i] *= -1;
  }
  return lowerBounds;
}

void interfaceConstraint::getBounds(patchSide ps, gsVector<> &bounds, index_t fixInfo){
    // Get the indicies corresponding the the patch.
    gsVector<index_t> I = getPatchSideIndicies(ps);

    // getTheControlPoints
    gsVector<index_t> J(I.size());
    getLocalPatchSideIndicies(ps,J);
    gsMatrix<> coefs = mp->patch(ps.patch).coefs();

    if (fixInfo != 0 and fixInfo != 1 and fixInfo != 2){
      GISMO_ERROR("In getBounds, fixInfo (k) has to be 0, 1 or 2");
    }

    // For each contol point do
    for(index_t k = 0; k < J.size(); k++){
      index_t ix = I[k];
      index_t iy = I[k]+n_controlpoints;
      if (fixInfo != 2){ // Fix x when k (fixInfo) is 0 or 1.
        bounds[ix] = mp->patch(ps.patch).coef(J[k],0);
      }
      if (fixInfo != 1){ // Fix y when k (fixInfo) is 0 or 2.
        bounds[iy] = mp->patch(ps.patch).coef(J[k],1);
      }
    }
}

void interfaceConstraint::freeBoundary(index_t p, index_t bnd, index_t k){
  // Method to free components of control points at boundaries (e.g. x coordinates)
  // for patch p, boundary::side bnd and freeing k
  // k = 0 : all fixed, k = 1 : x fixed, k = 2 : yfixed

  GISMO_ENSURE(k > -1 and k < 3,"Error in freeBoundary, k has to be 0, 1 or 2");

  index_t i = findBoundaryIndex(p,bnd);
  bndFixInfo[i] = k;
}

index_t interfaceConstraint::findBoundaryIndex(index_t p, index_t bnd){
  for( index_t i = 0; i < mp->nBoundary(); i++){
    patchSide ps = mp->boundaries()[i];
    // gsInfo << "\n\n" << ps.patch << ", " << ps.index() << std::flush;
    if (ps.patch == p and ps.index() == bnd){
      return i;
    }
  }
  GISMO_ERROR("In findBoundaryIndex: side " + std::to_string(bnd) + " on patch " + std::to_string(p) + " is not a Boundary for multipatch... ");
}
