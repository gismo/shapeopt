#ifndef GSHARMONIC_H
#define GSHARMONIC_H
using namespace gismo;

#include "gsHarmonic.h"
#include "gsOptParamMethod.h"

class gsHarmonic: public gsOptParamMethod {
public:
    // Use constructors from gsOptParamMethod, needs C++11 or newer for this..
    using gsOptParamMethod::gsOptParamMethod;

    // evaluation of objective
    real_t evalObj() const;

    gsVector<> gradObj() const;

    gsMatrix<> hessObj(gsMatrix<> &hessObjTagged) const;

    void setLambdas(real_t l1, real_t l2){ lambda_1 = l1; lambda_2 = l2;}

public:
    real_t lambda_1 = 1;
    real_t lambda_2 = 1;

};

#endif //GSHARMONIC_H
