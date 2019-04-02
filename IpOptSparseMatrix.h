#ifndef IPOPTSPARSEMATRIX_H
#define IPOPTSPARSEMATRIX_H
using namespace gismo;

class IpOptSparseMatrix{
public:
  void concatenate(IpOptSparseMatrix m2, std::string str);
  IpOptSparseMatrix(gsSparseMatrix<> mat);
  IpOptSparseMatrix(gsMatrix<> &mat, real_t tol=1e-10);
  IpOptSparseMatrix(){};
  const std::vector<index_t> & rows() const { return m_rows; }
  const std::vector<index_t> & cols() const { return m_cols; }
  const index_t & nnz() const { return m_nnz; }
  const index_t & nrows() const { return m_nrows; }
  const index_t & ncols() const { return m_ncols; }
  const gsVector<> & values() const { return m_values; }
  gsMatrix<> asDense();
  void printAsBlock();
  void print();
  void printSize();
  void swap(IpOptSparseMatrix smat);

  IpOptSparseMatrix& operator-(){ m_values *= -1; return *this; }

private:
  index_t m_nnz;
  index_t m_nrows;
  index_t m_ncols;
  gsVector<> m_values;
  std::vector<index_t> m_rows;
  std::vector<index_t> m_cols;

};


#endif //IPOPTSPARSEMATRIX_H
