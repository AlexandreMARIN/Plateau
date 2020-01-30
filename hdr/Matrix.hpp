#ifndef MATRIX_HPP
#define MATRIX_HPP



#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <list>

/* for displaying real vectors, represented by std::vector<double> */
std::ostream& operator<<(std::ostream&, const std::vector<double>&);


/**************Classes "Matrix"******************/


class Matrix{
protected:

  int nr, nc;

public:

  Matrix() = default;
  Matrix(const Matrix&) = default;
  Matrix(Matrix&&) = default;
  Matrix(int, int);
  virtual ~Matrix() = default;

  virtual Matrix& operator=(const Matrix&) = default;
  virtual Matrix& operator=(Matrix&&) = default;


  int get_nr() const;
  int get_nc() const;

  /*product of a matrix by a vector*/
  virtual void MvProd(const std::vector<double>&, std::vector<double>&) const = 0;

  virtual Matrix* clone() const = 0;

};

/********************************/
/*
  Class AtomicMatrix

  These matrices have entries
  we can access/modify.
*/
/********************************/

class AtomicMatrix : public Matrix{

public:

  AtomicMatrix() = default;
  AtomicMatrix(const AtomicMatrix&) = default;
  AtomicMatrix(AtomicMatrix&&) = default;
  AtomicMatrix(int, int);
  virtual ~AtomicMatrix() = default;

  virtual AtomicMatrix& operator=(const AtomicMatrix&) = default;
  virtual AtomicMatrix& operator=(AtomicMatrix&&) = default;

  /*for reading entries*/
  virtual double operator()(int, int) const = 0;
  virtual double& operator()(int, int) = 0;


  virtual void MvProd(const std::vector<double>&, std::vector<double>&) const = 0;

  void cg(std::vector<double>&, const std::vector<double>&) const;//conjugate gradient

  virtual AtomicMatrix* clone() const = 0;

};


/*************** Class DenseMatrix ******************/

class DenseMatrix : public AtomicMatrix{

  double* val;

public:

  DenseMatrix() = delete;
  DenseMatrix(int, int);
  DenseMatrix(const DenseMatrix&);
  DenseMatrix(DenseMatrix&&);
  ~DenseMatrix();

  DenseMatrix& operator=(const DenseMatrix&);
  DenseMatrix& operator=(DenseMatrix&&);
  double operator()(int, int) const;
  double& operator()(int, int);
  DenseMatrix operator()(const std::vector<int>&, const std::vector<int>&) const;

  /*diag(alpha, d) sets a diagonal of the matrix*/
  /*
    -if d=0, the diagonal will be filled with alpha
    -if d>0, the d-th diagonal under the diagonal will be filled with alpha
    -...
  */
  void resize(int, int);
  void diag(double, int = 0);
  void LoadFromFile(const std::string&);


  void MvProd(const std::vector<double>&, std::vector<double>&) const;
  void LUSolve(std::vector<double>&, const std::vector<double>&) const;


  DenseMatrix* clone() const override;

  friend std::ostream& operator<<(std::ostream&, const DenseMatrix&);

};


/*************** Class CSR ******************/

class CSR : public AtomicMatrix{

  int nnz;//number of non-zeros
  std::vector<std::list<std::pair<int, double> >::iterator> row;
  std::list<std::pair<int, double> > colval;

public:

  CSR() = default;
  CSR(int, int);
  CSR(const CSR&);
  CSR(CSR&&) = default;
  ~CSR() = default;


  CSR& operator=(const CSR&);
  CSR& operator=(CSR&&) = default;
  double operator()(int, int) const;
  double& operator()(int, int);

  void MvProd(const std::vector<double>&, std::vector<double>&) const;

  CSR* clone() const override;

  friend std::ostream& operator<<(std::ostream&, const CSR&);

};



#endif
