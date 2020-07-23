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

    virtual real_t evalObj(index_t p) const;

    // get gradient wrt. all control points
    virtual gsVector<> gradAll(gsDofMapper &space_mapper) const;

    // Returns hessian wrt. all control points
    virtual gsMatrix<> hessAll(gsDofMapper &space_mapper) const;

    void computeWinslowPerPatch();

    void addCorners(){ m_addCorners = true; };
    void setAlpha(real_t a ){ m_alpha = a; };

public:

	typedef memory::unique_ptr<gsWinslowPow> uPtr;
	typedef memory::shared_ptr<gsWinslowPow> Ptr;

public:

    real_t m_alpha = 0.01;
    bool m_addCorners = false;

};

#endif //GSWINSLOWPOW_H
