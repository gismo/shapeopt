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
    // Default empty constructor
    gsIpOptSparseMatrix(){};

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

    // Returns the matrix as dense
    gsMatrix<> asDense();

    // Accessors
    const std::vector<index_t> & rows() const { return m_rows; }
    const std::vector<index_t> & cols() const { return m_cols; }
    const index_t & nnz() const { return m_nnz; }
    const index_t & nrows() const { return m_nrows; }
    const index_t & ncols() const { return m_ncols; }
    const gsVector<> & values() const { return m_values; }

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
