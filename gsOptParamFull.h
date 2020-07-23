#ifndef GSOPTPARAMFULL_H
#define GSOPTPARAMFULL_H
using namespace gismo;

#include "gsOptParam.h"

class gsOptParamFull: public gsOptParam {
public:

    // use constructors from gsOptParam
    using gsOptParam::gsOptParam;

    gsVector<> getGradAllFromInterface(boundaryInterface interface) const;


    // Evaluation of the objective, using the design contained in m_mp
    real_t evalObj() const ;

    // Evaluation of the gradient wrt all cps, needed for gsShapeOptWithReg
    gsVector<> gradAll() const ;

    void setupMappers();

public:

	typedef memory::shared_ptr<gsOptParamFull> Ptr;
	typedef memory::unique_ptr<gsOptParamFull> uPtr;

public:

    std::vector< gsDofMapper > m_mappers_old;
    gsVector<index_t> patchShift;

};



# endif //GSOPTPARAMFULL_H
