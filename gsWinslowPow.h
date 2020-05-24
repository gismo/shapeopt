#ifndef GSWINSLOWPOW_H
#define GSWINSLOWPOW_H
using namespace gismo;

#include "gsWinslow.h"

class gsWinslowPow: public gsWinslow {
public:
    // Use constructors from gsWinslow, needs C++11 or newer for this..
    using gsWinslow::gsWinslow;

    // evaluation of objective
    virtual real_t evalObj() const;

    // get gradient wrt. all control points
    virtual gsVector<> gradAll(gsDofMapper &space_mapper) const;

    // Returns hessian wrt. all control points
    virtual gsMatrix<> hessAll(gsDofMapper &space_mapper) const;

public:

	typedef memory::unique_ptr<gsWinslowPow> uPtr;
	typedef memory::shared_ptr<gsWinslowPow> Ptr;

};

#endif //GSWINSLOWPOW_H
