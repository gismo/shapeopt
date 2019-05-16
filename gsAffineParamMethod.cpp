#include <gismo.h>
#include "gsParamMethod.h"
#include "gsAffineParamMethod.h"
using namespace gismo;

gsAffineParamMethod::gsAffineParamMethod(gsMultiPatch<>* mpin): gsParamMethod(mpin)
{
};

gsAffineParamMethod::gsAffineParamMethod(gsMultiPatch<>* mpin, std::vector< gsDofMapper > mappers):
        gsParamMethod(mpin, mappers)
{
    // computeMap();
};

// We linearize around 0 to obtain an affine map. The method doesnt make sense at 0,
// since all cps of patch will be 0.
// But it doesnt matter since we dont want to find good parametrizations at 0
// but just map out the affine function...
void gsAffineParamMethod::computeMap()
{
    // Compute b
    gsVector<> zers;
    zers.setZero(n_tagged);
    b = getUpdate(zers);

    // Compute A
    A.setZero(n_free,n_tagged);
    for(index_t i = 0; i < n_tagged; i++){
        gsVector<> ei;
        ei.setZero(n_tagged);
        ei[i] = 1;
        gsVector<> ci = getUpdate(ei) - b; // Substract value at zero since it is an affine function

        for(index_t j = 0; j < n_free; j++){
            A(j,i) = ci[j];
        }
    }
}

void gsAffineParamMethod::update(gsVector<> x)
{
    gsVector<> free_cps = b + A*x;
    // Update free and tagge DoFs
    updateFreeAndTagged(free_cps,x);
}

void gsAffineParamMethod::update()
{
    gsVector<> free_cps = getUpdate(getTagged());
    // Update free DoFs
    updateFree(free_cps);
}

// FIXIT.   Compute map has to be called before this method..
//          can we fix this? or at least
gsMatrix<> gsAffineParamMethod::jacobUpdate(gsVector<> x)
{
    // result doesnt depend on x..
    return A;
}
