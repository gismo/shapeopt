#ifndef GSMODLIAO_H
#define GSMODLIAO_H
using namespace gismo;

#include "gsOptParamMethod.h"

class gsModLiao: public gsOptParamMethod {
public:
    // Use constructors from gsOptParamMethod, needs C++11
    using gsOptParamMethod::gsOptParamMethod;

    // evaluation of objective
    real_t evalObj() const;

    gsVector<> gradObj() const;

    gsMatrix<> hessAll(gsDofMapper &space_mapper) const;

public:

	typedef memory::unique_ptr<gsModLiao> uPtr;
	typedef memory::shared_ptr<gsModLiao> Ptr;

};

#endif //GSMODLIAO_H
