#ifndef GSLIAO_H
#define GSLIAO_H
using namespace gismo;

#include "gsOptParamMethod.h"

class gsLiao: public gsOptParamMethod {
public:
    // Use constructors from gsOptParamMethod, needs C++11 or newer for this..
    using gsOptParamMethod::gsOptParamMethod;

    real_t evalObj() const;

    gsVector<> gradObj() const;

    gsMatrix<> hessAll(gsDofMapper &space_mapper) const;

};

#endif //GSLIAO_H
