#ifndef GSOPTPARAMFULL_H
#define GSOPTPARAMFULL_H
using namespace gismo;

#include "gsOptParam.h"

class gsOptParamFull: public gsOptParam {
public:

    // use constructors from gsOptParam
    using gsOptParam::gsOptParam;

    // Evaluation of the objective, using the design contained in m_mp
    real_t evalObj() const ;

    // Evaluation of the gradient wrt all cps, needed for gsShapeOptWithReg
    gsVector<> gradAll() const ;

    void setupMappers();

public:

	typedef memory::shared_ptr<gsOptParamFull> Ptr;
	typedef memory::unique_ptr<gsOptParamFull> uPtr;

public:

};



# endif //GSOPTPARAMFULL_H
