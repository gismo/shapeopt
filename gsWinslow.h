#ifndef GSWINSLOW_H
#define GSWINSLOW_H
using namespace gismo;

#include "gsOptParamMethod.h"

class gsWinslow: public gsOptParamMethod {
public:
    // Use constructors from gsOptParamMethod, needs C++11 or newer for this..
    using gsOptParamMethod::gsOptParamMethod;

    // Constructor with input flag that tells whether to check for infinite value of objective
    // \a checkForInf_eps will denote the value of detJ where objective will be set to 0
    gsWinslow(memory::shared_ptr<gsMultiPatch<>> mpin, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps);

    gsWinslow(memory::shared_ptr<gsMultiPatch<>> mpin, std::vector< gsDofMapper > mappers, bool use_dJC, bool useTensorStructureforDJC, bool checkForInf, real_t checkForInf_eps);

    // evaluation of objective
    virtual real_t evalObj() const;

    gsVector<> gradObj() const;

    // Returns hessian wrt. all variables
    gsMatrix<> hessAll() const;
    virtual gsMatrix<> hessAll(gsDofMapper &space_mapper) const;

    gsVector<> gradAll() const;
    virtual gsVector<> gradAll(gsDofMapper &space_mapper) const;

    gsVector<> gradObj(gsVector<> &gradObjTagged) const;

public:

	typedef memory::unique_ptr<gsWinslow> uPtr;
	typedef memory::shared_ptr<gsWinslow> Ptr;

public:
    // Flags that determines whether objective should be set to inf if detJ is smaller that m_checkForInf_eps in a gauss point.
    bool m_checkForInf = false;
    real_t m_checkForInf_eps = 0;

    real_t lambda = 1;
};

#endif //GSWINSLOW_H
