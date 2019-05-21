/** @file gsIpOptSparseMatrix.h

@brief Implements a sparse matrix in the format expected by IpOpt method jacobCon_into

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Limkilde, A. Mantzaflaris
*/
#ifndef GSIPOPTSPARSEMATRIX_H
#define GSIPOPTSPARSEMATRIX_H
using namespace gismo;

class gsIpOptSparseMatrix
{
public:
    // default empty constructer
    gsIpOptSparseMatrix(){};

    // Constructs empty (nxm) matrix
    gsIpOptSparseMatrix(index_t n, index_t m);

    // Constructs from a gsSparseMatrix using its sparsity pattern
    gsIpOptSparseMatrix(gsSparseMatrix<> mat);

    // Constructs from a gsMatrix, using sparsity pattern defined by the tolerance.
    // Setting tolorance < 0 (eg. -1) returns a dense matrix
    //FIXIT change default tolerance to -1 (dense), and check all other code is not ruined
    gsIpOptSparseMatrix(gsMatrix<> &mat, real_t tol=1e-10);

    // Operations with other gsIpOptSparseMatrices
    void concatenate(gsIpOptSparseMatrix m2, std::string str);
    void swap(gsIpOptSparseMatrix smat);
    gsIpOptSparseMatrix& operator-(){ m_values *= -1; return *this; }

    // FIXIT: Expensive!!
    // gsIpOptSparseMatrix& operator+=(gsIpOptSparseMatrix const& mat)
    // {
    //     std::vector<real_t> v;
    //     m_rows.clear();
    //     m_cols.clear();
    //
    //     for(index i = 0; i < m_nrows(); i++){
    //         for(index_t j = 0; j < m_ncols; j++){
    //             for(index_t k = 0; k < )
    //         }
    //     }
    //     // For each non zero value in this
    //     for(index_t j = 0; j < mat.nnz(); i++){
    //         for(index_t i = 0; i < m_nnz; i++){
    //             // If both has a value
    //             if (m_rows[i] = mat.rows(j) && m_cols[i] == mat.cols(j)){
    //                 m_values[i] += mat.values(i);
    //             } else {
    //                 m_rows.push_back(mat.rows(j));
    //                 m_cols.push_back(mat.cols(j));
    //                 m_nnz += 1;
    //
    //                 gsVector<> newVals(m_nnz);
    //                 newVals << m_values, mat.values(j);
    //             }
    //         }
    //     }
    //     return *this;
    // }

    // Returns the matrix as dense
    gsMatrix<> asDense();

    // Returns the matrix as gsSparseMatrix
    gsSparseMatrix<> asSparse();

    // Accessors
    // FIXIT: think about the naming such that it is similar to the gsMatrx
    std::vector<index_t> & rows() { return m_rows; }
    std::vector<index_t> & cols() { return m_cols; }
    const index_t & nnz() const { return m_nnz; }
    index_t & nrows() { return m_nrows; }
    index_t & ncols() { return m_ncols; }
    const gsVector<> & values() const { return m_values; }

    const index_t rows(index_t i) const { return m_rows[i]; }
    const index_t cols(index_t i) const { return m_cols[i]; }
    const real_t values(index_t i) const { return m_values[i]; }

    // Setters
    void setValues(gsVector<> values) { m_values = values; }
    void setNnz(index_t nnz) { m_nnz = nnz; }

    // Extract a column of the matrix
    gsIpOptSparseMatrix col(index_t i);

    // Print
    void printAsBlock();
    void print();
    void printSize();

private:
    index_t m_nnz;
    index_t m_nrows;
    index_t m_ncols;
    gsVector<> m_values;
    std::vector<index_t> m_rows;
    std::vector<index_t> m_cols;

};


#endif //IPOPTSPARSEMATRIX_H
