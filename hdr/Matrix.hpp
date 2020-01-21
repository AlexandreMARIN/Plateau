#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>

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

  /**************************************/
  /*To compare two sorts of conjugate****/
  /*gradients, we use them on the matrix*/
  /*and save norms of residual vectors***/
  /*in a file in Json format.************/
  /**************************************/
  std::string resvec_to_json(double, int) const;

  virtual void MvProd(const std::vector<double>&, std::vector<double>&) const = 0;
  virtual void LUSolve(std::vector<double>&, const std::vector<double>&) const = 0;
  void cg(std::vector<double>&, const std::vector<double>&) const;
  void pcg(std::vector<double>&, const std::vector<double>&, double) const;
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
  int* row, *col;
  double* val;

public:

  CSR() = default;
  CSR(int, int);
  CSR(const DenseMatrix&);
  CSR(const CSR&);
  CSR(CSR&&);
  ~CSR();

  operator DenseMatrix() const;

  CSR& operator=(const CSR&);
  CSR& operator=(CSR&&);
  double operator()(int, int) const;
  double& operator()(int, int);
  CSR operator()(const std::vector<int>&, const std::vector<int>&) const;


  void MvProd(const std::vector<double>&, std::vector<double>&) const;
  void LUSolve(std::vector<double>&, const std::vector<double>&) const;

  CSR* clone() const override;

  friend std::ostream& operator<<(std::ostream&, const CSR&);

  /*these methods quickly build sparse matrices*/
  static CSR laplacian(int);
  static CSR H(int);
};



#endif
