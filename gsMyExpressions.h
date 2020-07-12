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

template <typename E1, typename E2>
class my_collapse_expr : public _expr<my_collapse_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = 0, ColBlocks = 0};
    enum { Space = E1::Space || E2::Space};

    typedef typename E1::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    my_collapse_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) { }

    //EIGEN_STRONG_INLINE MatExprType
    const gsMatrix<Scalar> &
    eval(const index_t k) const
    {
        const index_t nb = rows();
        const MatExprType tmpA = _u.eval(k);
        const MatExprType tmpB = _v.eval(k);

        /*
        gsDebugVar(nb);
        gsDebugVar(tmpA.rows());
        gsDebugVar(tmpA.cols());
        gsDebugVar(tmpB.rows());
        gsDebugVar(tmpB.cols());
        gsDebugVar(_v.rows());
        gsDebugVar(E1::ColBlocks);
        gsDebugVar(E2::ColBlocks);
        */

        if (E1::ColBlocks)
        {
            const index_t ur = _v.rows();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).transpose().noalias() = tmpA.middleCols(i*ur,ur) * tmpB;
            }
        }
        else if (E2::ColBlocks)
        {
            const index_t ur = _u.cols();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).noalias() = tmpA * tmpB.middleCols(i*ur,ur);
            }
        }

        return res;
    }

    index_t rows() const { return E1::ColBlocks ? _u.cols() / _v.rows() : _v.cols() / _u.cols() ; }
    index_t cols() const { return E1::ColBlocks ? _v.rows()  : _u.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    static constexpr bool rowSpan() { return E1::Space ? E1::rowSpan() : E2::rowSpan(); }
    static bool colSpan() { return E2::Space ? E2::colSpan() : E1::colSpan(); }

    const gsFeSpace<Scalar> & rowVar() const
    { return E1::ColBlocks ? _u.rowVar() : _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {
        GISMO_ERROR("none");
    }

    void print(std::ostream &os) const { _u.print(os); os<<"my~"; _v.print(os); }
};

// Multi-matrix collapsed by a vector
template <typename E1, typename E2> //EIGEN_STRONG_INLINE
//collapse_expr<E1,E2> const  operator&(<E1> const& u, _expr<E2> const& v)
my_collapse_expr<E1,E2> my_collapse( _expr<E1> const& u, _expr<E2> const& v)
{ return my_collapse_expr<E1, E2>(u, v); }

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

/*
   Expression for the derivate of the normal vector - Local to ASGL
 */
template<class E>
class nvDeriv_expr : public _expr<nvDeriv_expr<E> >
{
    typename E::Nested_t _u;
    gsGeometryMap<typename E::Scalar> _G;

    public:
    enum {Space = E::Space };

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    nvDeriv_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

    // FIXIT: Maybe this stuff can be more efficient. Look again, once validated.
    MatExprType eval(const index_t k) const
    {
        const index_t dd = _u.source().domainDim();
        index_t n = _u.rows();
        res.setZero(n, dd);

        gsMatrix<> jac = _u.data().values[1].col(k).transpose()
            .blockDiag(_u.dim()); 

        // FIXIT make dimension independent! Maybe ask Angelos
        const gsAsConstMatrix<Scalar, 3, 3> jacG(_G.data().values[1].col(k).data(), 3, 3);


        real_t sgn = sideOrientation( _G.data().side );
        index_t dir = _G.data().side.direction();

        
        //gsDebugVar(dd);
        //gsDebugVar(n);
        //gsDebugVar(_u.cols());
        //gsDebugVar(jac.rows());
        //gsDebugVar(jac.cols());
        for ( index_t j = 0; j < n; j++)
        {
            gsMatrix<> _jacT = jac.block(0,dd*j,dd,dd).transpose();
            const gsAsConstMatrix<Scalar, 3, 3> jacT_j( _jacT.data(), 3, 3 );

            index_t alt_sgn = sgn * ( jacG.determinant()<0 ? -1 : 1);

            typename gsMatrix<Scalar, 3, 3>::FirstMinorMatrixType minor, minorG;
            for ( index_t i = 0; i < dd; i++)
            {

                jacT_j.firstMinor(dir, i, minor);
                jacG.firstMinor(dir, i, minorG);

                gsMatrix<Scalar, 2, 2> tmpG(minorG);

                res(j,i) = alt_sgn*(tmpG.adjugate()*minor).trace();// Derivative of det.!
                alt_sgn *= -1;
            }
        }


        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t rows() const { 
        return _u.data().values[1].rows();
    }

    index_t cols() const { return _u.source().domainDim(); }

    static constexpr bool rowSpan() {return true; }
    static bool colSpan() {return false;}

    void setFlag() const
    {
        _u.data().flags |= NEED_DERIV;
        _G.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(& _u.source());
        evList.push_sorted_unique(& _G.source());

        _u.data().flags |= NEED_DERIV;
        _G.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void print(std::ostream &os) const { os << "nvDeriv("; _u.print(os); os <<")"; }
};

template<class T>
class mynormal_expr : public _expr<mynormal_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    mynormal_expr(const gsGeometryMap<T> & G) : _G(G) { }

    MatExprType eval(const index_t k) const
    {
        gsVector<> Exact_Val = _G.data().outNormals.col(k);

        const index_t dd = 3;

        // FIXIT make dimension independent! Maybe ask Angelos
        const gsAsConstMatrix<Scalar, 3, 3> jacG(_G.data().values[1].col(k).data(), 3, 3);

        real_t sgn = sideOrientation( _G.data().side );
        index_t dir = _G.data().side.direction();
        
        index_t alt_sgn = sgn * ( jacG.determinant()<0 ? -1 : 1);

        gsVector<> out;
        out.setZero(3);

        typename gsMatrix<Scalar, 3, 3>::FirstMinorMatrixType minor, minorG;
        for ( index_t i = 0; i < dd; i++)
        {
            jacG.firstMinor(dir, i, minorG);

            gsMatrix<Scalar, 2, 2> tmpG(minorG);

            out[i] = alt_sgn*tmpG.determinant();// Derivative of det.!
            alt_sgn *= -1;
        }


        return out;
    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    static constexpr bool rowSpan() {return false;}
    static bool colSpan() {return false;}

    void setFlag() const { _G.data().flags |= NEED_OUTER_NORMAL; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    // Normalized to unit length
    normalized_expr<mynormal_expr<T> > normalized()
    { return normalized_expr<mynormal_expr<T> >(*this); }

    void print(std::ostream &os) const { os << "myNv("; _G.print(os); os <<")"; }
};

/*
   Expression for the derivate of the normal vector - Local to ASGL
 */
template<class E>
class collapsedJac_expr : public _expr<collapsedJac_expr<E> >
{
    typename E::Nested_t _u;

    public:
    enum {Space = E::Space };

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    collapsedJac_expr(const E & u) : _u(u) { }

    MatExprType eval(const index_t k) const
    {
        const index_t dd = _u.source().domainDim();
        index_t n = _u.rows();
        res.setZero(n, dd);

        gsMatrix<> jac = _u.data().values[1].col(k).transpose()
            .blockDiag(_u.dim()); 

        //gsDebugVar(dd);
        //gsDebugVar(n);
        //gsDebugVar(_u.cols());
        //gsDebugVar(jac.rows());
        //gsDebugVar(jac.cols());
        for ( index_t j = 0; j < n; j++)
        {
            res.row(j) = jac.block(0,dd*j,dd,dd).diagonal();
        }

        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t rows() const { 
        return _u.data().values[1].rows() / cols();
    }

    index_t cols() const { return _u.source().domainDim(); }

    static constexpr bool rowSpan() {return true; }
    static bool colSpan() {return false;}

    void setFlag() const
    {
        _u.data().flags |= NEED_DERIV;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(& _u.source());

        _u.data().flags |= NEED_DERIV;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void print(std::ostream &os) const { os << "collapsedJac("; _u.print(os); os <<")"; }
};

/// The derivative of the jacobian of a geometry map with respect to a coordinate.
template<class E> EIGEN_STRONG_INLINE
dJacdc_expr<E> dJacdc(const E & u, index_t c) { return dJacdc_expr<E>(u,c); }

/// The transpose of Jacobian matrix of a FE variable
template<class E> EIGEN_STRONG_INLINE
jacT_expr<E> jacT(const E & u) { return jacT_expr<E>(u); }

/// The (outer pointing) boundary normal of a geometry map
template<class T> EIGEN_STRONG_INLINE
mynormal_expr<T> myNv(const gsGeometryMap<T> & u) { return mynormal_expr<T>(u); }

/// The derivate of the normal vector - Local to ASGL
template<class E> EIGEN_STRONG_INLINE
nvDeriv_expr<E> nvDeriv(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return nvDeriv_expr<E>(u,G); }


/// The derivate of the normal vector - Local to ASGL
template<class E> EIGEN_STRONG_INLINE
collapsedJac_expr<E> collapsedJac(const E & u) { return collapsedJac_expr<E>(u); }


} // namespace expr

} //namespace gismo
