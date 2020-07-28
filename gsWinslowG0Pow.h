#ifndef GSWINSLOWG0POW_H
#define GSWINSLOWG0POW_H
using namespace gismo;

#include "gsWinslow.h"

class gsWinslowG0Pow: public gsWinslow {
public:
    // Use constructors from gsWinslow, needs C++11 or newer for this..
    using gsWinslow::gsWinslow;

    // evaluation of objective
    virtual real_t evalObj() const;

    // get gradient wrt. all control points
    virtual gsVector<> gradAll(gsDofMapper &space_mapper) const;

    void addCorners(){ m_addCorners = true; };
    void setAlpha(real_t a ){ m_alpha = a; };

    void setMp0(gsMultiPatch<> mp){ m_mp0 = mp; };

public:

	typedef memory::unique_ptr<gsWinslowG0Pow> uPtr;
	typedef memory::shared_ptr<gsWinslowG0Pow> Ptr;

public:

    gsMultiPatch<> m_mp0;

    real_t m_alpha = 0.01;
    bool m_addCorners = false;

};

#endif //GSWINSLOWG0POW_H
