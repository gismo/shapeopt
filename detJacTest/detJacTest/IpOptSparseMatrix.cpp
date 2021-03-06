#include <gismo.h>
#include "IpOptSparseMatrix.h"
using namespace gismo;

IpOptSparseMatrix::IpOptSparseMatrix(gsSparseMatrix<> mat){
  m_nnz = mat.nonZeros();
  m_nrows = mat.rows();
  m_ncols = mat.cols();
  m_values.setZero(m_nnz);
  index_t i = 0;
  for (int k=0; k<mat.outerSize(); ++k)
  for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
  {
    m_values[i++] = it.value();
    m_rows.push_back(it.row());   // row index
    m_cols.push_back(it.col());   // col index
  }
}

IpOptSparseMatrix::IpOptSparseMatrix(gsMatrix<> &mat,real_t tol){
  m_nnz = 0;

  auto zmat = mat.cwiseAbs().array() > tol;
  m_nnz = zmat.count();

  m_nrows = mat.rows();
  m_ncols = mat.cols();
  m_values.setZero(m_nnz);
  index_t i = 0;
  for(index_t k = 0; k < m_nrows; k++){
    for(index_t l = 0; l < m_ncols; l++){
      if (zmat(k,l)){
        m_values[i] = mat(k,l);
        m_rows.push_back(k);
        m_cols.push_back(l);
        i++;
      }
    }
  }
}

// Method to concatenate two matricies A,B, in the fashion
//  str=row : [A B]
//  str=col : [A
//             B]
//  str=diag: [A 0
//             0 B]
void IpOptSparseMatrix::concatenate(IpOptSparseMatrix m2, std::string str = "diag" ){
  index_t nnz = m_nnz + m2.nnz();
  gsVector<> values;
  values.setZero(nnz);

  index_t addToRows = 0, addToCols = 0;
  if (str.compare("diag") == 0) {
      addToRows += m_nrows;
      addToCols += m_ncols;
      m_nrows += m2.nrows();
      m_ncols += m2.ncols();
  } else if (str.compare("row") == 0) {
    addToCols += m_ncols;
    m_ncols += m2.ncols();
  } else if (str.compare("col") == 0) {
    addToRows += m_nrows;
    m_nrows += m2.nrows();
  } else {
    GISMO_ERROR("IpOptSparseMatrix concatenate constructer expects str to be equal 'diag', 'row' or 'col'");
  }

  index_t i = 0;
  for(index_t k = 0; k < m_nnz; k++){
    values[i] = m_values[k];
    i++;
  }
  for(index_t k = 0; k < m2.nnz(); k++){
    values[i] = m2.values()[k];
    m_rows.push_back(m2.rows()[k] + addToRows);
    m_cols.push_back(m2.cols()[k] + addToCols);
    i++;
  }
  m_nnz = nnz;
  m_values = values;
}

void IpOptSparseMatrix::print(){
  for(index_t i = 0; i < m_nnz; i++){
    gsInfo << "(" << m_rows[i] << ", " << m_cols[i] << "): " << m_values[i] << "\n";
  }
}

void IpOptSparseMatrix::printSize(){
  gsInfo << "(" << m_nrows << ", " << m_ncols << ")\n";
}

void IpOptSparseMatrix::printAsBlock(){
  gsInfo << this->asDense() << "\n";
}

gsMatrix<> IpOptSparseMatrix::asDense(){
  gsMatrix<> out;
  out.setZero(m_nrows,m_ncols);
  for(index_t i = 0; i < m_nnz; i++){
    out(m_rows[i],m_cols[i]) = m_values[i];
  }
  return out;
}

void IpOptSparseMatrix::swap(IpOptSparseMatrix smat){
  this->m_nrows = smat.nrows();
  this->m_ncols = smat.ncols();
  this->m_nnz = smat.nnz();
  this->m_values = smat.values();
  this->m_rows = smat.rows();
  this->m_cols = smat.cols();
}
