#ifndef GSWINSLOW_H
#define GSWINSLOW_H
using namespace gismo;

#include "gsOptParamMethod.h"

class gsWinslow: public gsOptParamMethod {
public:
    // Use constructors from gsOptParamMethod, needs C++11 or newer for this..
    using gsOptParamMethod::gsOptParamMethod;

    // evaluation of objective
    real_t evalObj() const;

    gsVector<> gradObj() const;

    gsMatrix<> hessObj(gsMatrix<> &hessObjTagged) const;

};

#endif //GSWINSLOW_H
