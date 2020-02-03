/** @file gsMyExpressions.h

    @brief Defines different expressions that are not in the gsExpressions.h file already

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Limkilde
*/

#pragma once

#include <gismo.h>

namespace gismo
{

// Forward declaration in gismo namespace
template<class T> class gsExprHelper;

namespace expr
{

#if(__cplusplus >= 201402L) // c++14
#  define MatExprType  auto
#  define AutoReturn_t auto
//#elif(__cplusplus >= 201103L) // c++11
//note: in c++11 auto-return requires -> decltype(.)
#else // 199711L
#  define MatExprType typename gsMatrix<Scalar>::constRef
#  define AutoReturn_t typename util::conditional<ScalarValued,Scalar,MatExprType>::type
#endif

/*
   Expression for the derivative of the jacobian of a spline geometry map,
   with respect to the coordinate c.

   It returns a matrix with the gradient of u in row d.
 */
template<class E>
class dJacdc_expr : public _expr<dJacdc_expr<E> >
{
    typename E::Nested_t _u;

    public:
    enum{ Space = E::Space, ScalarValued = 0, ColBlocks = E::rowSpan()};

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    index_t _c;

    dJacdc_expr(const E & u, index_t c) : _u(u), _c(c)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t dd = _u.source().domainDim();
        index_t n = _u.rows();
        res.setZero(dd, dd*n);

        gsMatrix<Scalar> grad = _u.data().values[1].reshapeCol(k, dd, n);
        for(index_t i = 0; i < n; i++){
            res.row(_c).segment(i*dd,dd) = grad.col(i);
        }
        return res;
    }

    index_t rows() const
    {
        //return _u.data().values[0].rows();
        return _u.source().domainDim();
    }
    //index_t rows() const { return _u.data().actives.size(); }
    //index_t rows() const { return _u.rows(); }

    //index_t rows() const { return _u.source().targetDim() is wrong }
    index_t cols() const { return _u.source().domainDim()*_u.rows(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void print(std::ostream &os) const { os << "dJacdc("; _u.print(os); os <<")"; }
};

/*
   Expression for the transpose Jacobian matrix of a FE variable - Local to ASGL
 */
template<class E>
class jacT_expr : public _expr<jacT_expr<E> >
{
    typename E::Nested_t m_fev;

    public:
    enum {ColBlocks = E::rowSpan() };
    enum {Space = E::Space };

    typedef typename E::Scalar Scalar;

    jacT_expr(const E & _u)
    : m_fev(_u) { }

    MatExprType eval(const index_t k) const
    {
        // Dim x (numActive*Dim)
        // gsMatrix<> jacu = m_fev.data().values[1].col(k).transpose().blockDiag(m_fev.dim());
        // const index_t r = jacu.rows();
        // const index_t n = jacu.cols()/r;
        //
        // for(index_t i = 0; i < n; i++){
        //     jacu.middleCols(i*r,r).transposeInPlace();
        // }

        return m_fev.data().values[1].col(k).transpose().blockDiag(m_fev.dim()).blockTranspose(m_fev.dim());
    }

    const gsFeSpace<Scalar> & rowVar() const { return m_fev; }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t rows() const { return m_fev.dim(); }
    index_t cols() const
    {   //bug
        return m_fev.dim() * m_fev.data().actives.rows() * m_fev.data().dim.first;
    }

    static constexpr bool rowSpan() {return true; }
    static bool colSpan() {return false;}

    void setFlag() const
    {
        m_fev.data().flags |= NEED_DERIV;
        m_fev.data().flags |= NEED_ACTIVE;// rows() depend on this
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(& m_fev.source());
         m_fev.data().flags |= NEED_DERIV;
         m_fev.data().flags |= NEED_ACTIVE;// rows() depend on this
    }

    void print(std::ostream &os) const { os << "jacT("; m_fev.print(os);os <<")"; }
};

/// The derivative of the jacobian of a geometry map with respect to a coordinate.
template<class E> EIGEN_STRONG_INLINE
dJacdc_expr<E> dJacdc(const E & u, index_t c) { return dJacdc_expr<E>(u,c); }

/// The transpose of Jacobian matrix of a FE variable
template<class E> EIGEN_STRONG_INLINE
jacT_expr<E> jacT(const E & u) { return jacT_expr<E>(u); }


} // namespace expr

} //namespace gismo
