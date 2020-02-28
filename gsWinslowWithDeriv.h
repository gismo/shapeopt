#ifndef GSWINSLOWWITHDERIV_H
#define GSWINSLOWWITHDERIV_H
using namespace gismo;

#include "gsWinslow.h"
#include "gsAffineOptParamMethod.h"

class gsWinslowWithDeriv: public gsWinslow {
public:
    // Constructor, we only include the most general one here as of now
    gsWinslowWithDeriv(memory::shared_ptr<gsMultiPatch<>> mpin, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps);

    gsWinslowWithDeriv(memory::shared_ptr<gsMultiPatch<>> mpin, std::vector< gsDofMapper > mappers, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps);

    using gsOptParamMethod::update;

    bool checkForNegativeDetJ();

    real_t minDetJInGaussPts(index_t incPts = 0);

    // Method to update inner controlpoints, given tagged Dofs (x).
    // Checks that updating the tagged vertices will not lead to negative determinant
    // Returns false if the update failed
    bool update(gsVector<> x);

    // Method to get the jacobian of the update (free Dofs) with respect to x (the tagged Dofs)
    gsMatrix<> jacobUpdate(gsVector<> x);

public:

	typedef memory::shared_ptr<gsWinslowWithDeriv> Ptr;
	typedef memory::unique_ptr<gsWinslowWithDeriv> uPtr;

public:
    gsAffineOptParamMethod m_aff;

};

#endif //GSWINSLOWWITHDERIV_H
