/** @file gsAggFun.h

@brief  A function that can be used for constraint aggregation

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSAGGFUN_H
#define GSAGGFUN_H

using namespace gismo;

class gsAggFun: public gsFunction<real_t>{
public:

    // Constructs Identity function for targetDim() = domainDim() = \a n
    gsAggFun(index_t n, real_t alpha): m_n(n), m_alpha(alpha) {};

    // Evaluate function, now Identity
    void eval_into(const gsMatrix<>& u, gsMatrix<>& result) const {
        result.setZero(1,u.cols());
        for (index_t i = 0; i < u.cols(); i++){
            gsVector<> eax = (m_alpha*u.col(i)).array().exp();
            result(0,i) = eax.dot(u.col(i))/eax.sum();
        }
    };

    short_t domainDim() const { return m_n; };

    short_t targetDim() const { return 1; };


public:
    index_t m_n;
    real_t m_alpha;
};




#endif //GSAGGFUN_H
