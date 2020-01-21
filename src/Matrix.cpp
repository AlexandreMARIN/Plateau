#include "../include/Block.hpp"

#include <chrono>
#include <random>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>

using namespace std;
using namespace std::chrono;

/*
For displaying vectors of doubles
*/
ostream& operator<<(ostream& os, const vector<double>& v){
  os << "[ ";
  if(v.size()){
    os << setw(10) << v[0];
  }
  for(vector<double>::size_type i=1;i<v.size();i++){
    os << ", " << setw(10) << v[i];
  }
  os << " ]";

  return os;
}


/*****************************************/
/*************Class Matrix****************/
/*****************************************/

Matrix::Matrix(int r, int c) : nr(r), nc(c){}

int Matrix::get_nr() const{
  return nr;
}

int Matrix::get_nc() const{
  return nc;
}


/*****************************************/
/***********Class AtomicMatrix************/
/*****************************************/

AtomicMatrix::AtomicMatrix(int r, int c) : Matrix(r, c){}

void AtomicMatrix::cg(vector<double>& x, const vector<double>& b) const{
  assert(nc==nr && (vector<double>::size_type)nc==b.size());

  ConjGrad_.set_A(this);
  ConjGrad_.set_b(&b);
  ConjGrad_.solve();
  x = ConjGrad_.get_x();
}

void AtomicMatrix::pcg(vector<double>& x, const vector<double>& b, double epsilon) const{
  assert(nc==nr && (vector<double>::size_type)nc==b.size());

  ASM_.set_A(this);
  ASM_.set_b(&b);
  ASM_.set_tol(epsilon);
  ASM_.solve();
  x = ASM_.get_x();
}

string AtomicMatrix::resvec_to_json(double epsilon, int n_max) const{
  ostringstream ss;
  vector<double> x, b;
  high_resolution_clock::time_point tp = high_resolution_clock::now();
  linear_congruential_engine<long unsigned int, 16807ul, 0ul, 2147483647ul> lce{static_cast<linear_congruential_engine<long unsigned int, 16807ul, 0ul, 2147483647ul>::result_type>(tp.time_since_epoch().count())};
  uniform_real_distribution<> urd{};

  for(int i=0;i<nc;i++){
    b.push_back(urd(lce));
  }

  ConjGrad_.set_n_max(n_max);
  ConjGrad_.set_tol(epsilon);
  ASM_.set_n_max(n_max);
  cg(x, b);
  ss << "\n\t\t\"resvec_cg\" : ";
  ss << ConjGrad_.get_resvec();
  ss << ",\n\t\t\"resvec_pcg\" : ";
  pcg(x, b, epsilon);
  ss << ASM_.get_resvec();
  return ss.str();
}


/*****************************************/
/************Class DenseMatrix************/
/*****************************************/

DenseMatrix::DenseMatrix(int r, int c) : AtomicMatrix(r, c), val(new double[r*c]{}){}

DenseMatrix::DenseMatrix(const DenseMatrix& A) : AtomicMatrix(A.nr, A.nc), val(new double[nr*nc]){

  for(int i=0;i<nr*nc;i++){
    val[i] = A.val[i];
  }

}

DenseMatrix::DenseMatrix(DenseMatrix&& A) : AtomicMatrix(A.nr, A.nc), val(A.val){
  A.nr = 0;
  A.nc = 0;
  A.val = nullptr;
}

DenseMatrix::~DenseMatrix(){
  delete[] val;
}

DenseMatrix& DenseMatrix::operator=(const DenseMatrix& A){
  if(nr*nc!=A.nr*A.nc){
    delete[] val;
    val = new double[A.nr*A.nc];
  }

  nr = A.nr;
  nc = A.nc;

  for(int i=0;i<nr*nc;i++){
    val[i] = A.val[i];
  }

  return *this;
}

DenseMatrix& DenseMatrix::operator=(DenseMatrix&& A){
  delete[] val;
  nr = A.nr;
  nc = A.nc;
  val = A.val;
  A.nr = 0;
  A.nc = 0;
  A.val = nullptr;

  return *this;
}

double DenseMatrix::operator()(int i, int j) const{
  assert(i>=0 && j>=0 && i<nr && j<nc);
  return val[i*nc+j];
}

double& DenseMatrix::operator()(int i, int j){
  assert(i>=0 && j>=0 && i<nr && j<nc);
  return val[i*nc+j];
}

DenseMatrix DenseMatrix::operator()(const vector<int>& Ir, const vector<int>& Ic) const{
  DenseMatrix A(static_cast<int>(Ir.size()), static_cast<int>(Ic.size()));

  int i, j;
  for(i=0;(vector<int>::size_type)i<Ir.size();i++){
    for(j=0;(vector<int>::size_type)j<Ic.size();j++){
      A.val[i*static_cast<int>(Ic.size())+j] = val[Ir[i]*nc+Ic[j]];
    }
  }

  return A;
}


void DenseMatrix::diag(double alpha, int d){
  int aux = (nr<d+nc)?nr:d+nc;
  for(int i=(d>=0)?d:0;i<aux;i++){
    val[i*nc+i-d] = alpha;
  }
}

void DenseMatrix::LoadFromFile(const string& s){
  delete val;
  ifstream file(s);
  char buf[10];
  file.read(buf, 6);
  file >> nr;
  file >> nc;
  file.read(buf, 7);
  val = new double[nr*nc];
  for(int i=0;i<nr*nc;i++){
    file >> val[i];
  }
}

void DenseMatrix::MvProd(const vector<double>& x, vector<double>& b) const{

  int i, j;

  assert((vector<double>::size_type)nc==x.size());

  b.resize(x.size());

  for(i=0;i<nr;i++){
    b[i] = 0.0;
    for(j=0;j<nc;j++){
      b[i] += val[i*nc+j]*x[j];
    }
  }

}

void DenseMatrix::LUSolve(vector<double>& x, const vector<double>& b) const{

  assert(nc==nr && (vector<double>::size_type)nc==b.size());

  x.resize(b.size());

  int i, j, k;
  DenseMatrix LU(*this);

  //we compute the LU decomposition
  for(k=0;k<=nc-2;k++){
    for(i=1;i<nc-k;i++){
      LU.val[(k+i)*nc+k] /= LU.val[k*(nc+1)];
      for(j=1;j<nc-k;j++){
	LU.val[(k+i)*nc+k+j] -= LU.val[(k+i)*nc+k]*LU.val[k*nc+j+k];
      }
    }
  }

  //we solve LX = b
  x = b;
  for(i = 1;i<nc;i++){
    for(j = 0;j<i;j++){
      x[i] -= LU.val[i*nc+j]*x[j];
    }
  }

  //we solve Ux = X
  for(i = nc-1;i>=0;i--){
    for(j = i+1;j<nc;j++){
      x[i] -= LU.val[i*nc+j]*x[j];
    }
    x[i] /= LU.val[(nc+1)*i];
  }

}

DenseMatrix* DenseMatrix::clone() const{
  return new DenseMatrix(*this);
}


ostream& operator<<(ostream& os, const DenseMatrix& A){

  os << A.nr << " by " << A.nc << " DenseMatrix :\n";

  int i, j;
  for(i=0;i<A.nr;i++){
    os << "| ";
    for(j=0;j<A.nc;j++){
      os << setw(10) << A.val[i*A.nc+j] << " ";
    }
    os << "|\n";
  }

  os << "\n";

  return os;
}


/*****************************************/
/****************Class CSR****************/
/*****************************************/

CSR::CSR(int r, int c) : AtomicMatrix(r, c), nnz(0), row(new int[r]{}), col(nullptr), val(nullptr){}

CSR::CSR(const DenseMatrix& A) : AtomicMatrix(A), nnz(0), row(new int[nr]), col(nullptr), val(nullptr){

  int* col2 = new int[nr*nc];
  double* val2 = new double[nr*nc];
  int i, j;
  bool nullrow, nullprevrow;
  for(i=0;i<nr;i++){
    nullrow = true;
    for(j=0;j<nc;j++){
      if(A(i, j)!=0.0){
	nullrow = false;
	break;
      }
    }
    if(nullrow){
      if(i==0){
	row[i] = 0;
      }else{
	if(nullprevrow){
	  row[i] = row[i-1];
	}else{
	  row[i] = 1+nnz;
	}
      }
    }else{
      row[i] = nnz;
      col2[nnz] = j;
      val2[nnz] = A(i, j);
      for(j++, nnz++;j<nc;j++){
	if(A(i, j)!=0.0){
	  col2[nnz] = j;
	  val2[nnz++] = A(i, j);
	}
      }
    }
    nullprevrow = nullrow;
  }

  col = new int[nnz];
  val = new double[nnz];
  for(i=0;i<nnz;i++){
    col[i] = col2[i];
    val[i] = val2[i];
  }
  delete col2;
  delete val2;
}

CSR::CSR(const CSR& A) : AtomicMatrix(A), nnz(A.nnz), row(new int[A.nr]), col(new int[nnz]), val(new double[nnz]){
  int i;
  for(i=0;i<nr;i++){
    row[i] = A.row[i];
  }
  for(i=0;i<nnz;i++){
    col[i] = A.col[i];
    val[i] = A.val[i];
  }
}

CSR::CSR(CSR&& A) : AtomicMatrix(A), nnz(A.nnz), row(A.row), col(A.col), val(A.val){
  A.nnz = 0;
  A.nr = 0;
  A.nc = 0;
  A.row = nullptr;
  A.col = nullptr;
  A.val = nullptr;
}

CSR::~CSR(){
  delete row;
  delete col;
  delete val;
}

CSR::operator DenseMatrix() const{
  DenseMatrix A(nr, nc);
  int i, j;
  for(i=0;i+1<nr;i++){
    for(j=row[i];j<row[i+1];j++){
      A(i, col[j]) = val[j];
    }
  }
  //i==nr-1
  for(j=row[i];j<nnz;j++){
    A(i, col[j]) = val[j];
  }
  return A;
}

CSR& CSR::operator=(const CSR& A){
  if(nr!=A.nr){
    nr = A.nr;
    delete row;
    row = new int[nr];
  }

  if(nnz!=A.nnz){
    nnz = A.nnz;
    delete col;
    delete val;
    col = new int[nnz];
    val = new double[nnz];
  }

  int i;
  for(i=0;i<nr;i++){
    row[i] = A.row[i];
  }
  for(i=0;i<nnz;i++){
    col[i] = A.col[i];
    val[i] = A.val[i];
  }

  nc = A.nc;
  return *this;
}

CSR& CSR::operator=(CSR&& A){
  nnz = A.nnz;
  nr = A.nr;
  nc = A.nc;
  row = A.row;
  col = A.col;
  val = A.val;
  A.nnz = 0;
  A.nr = 0;
  A.nc = 0;
  A.row = nullptr;
  A.col = nullptr;
  A.val = nullptr;
  return *this;
}

double CSR::operator()(int i, int j) const{

  int k;

  if(i==nr-1){
    for(k=row[i];k<nnz;k++){
      if(j==col[k]){
	return val[k];
      }
      if(j<col[k]){
	return 0.0;
      }
    }
  }else{
    for(k=row[i];k<row[i+1];k++){
      if(j==col[k]){
	return val[k];
      }
      if(j<col[k]){
	return 0.0;
      }
    }
  }

  return 0.0;
}

double& CSR::operator()(int i, int j){

  int k;

  if(i==nr-1){
    for(k=row[i];k<nnz;k++){
      if(j==col[k]){
	return val[k];
      }
      if(j<col[k]){
	throw(invalid_argument("double& CSR::operator()(int, int) : bad indices\n"));
      }
    }
  }

  for(k=row[i];k<row[i+1];k++){
    if(j==col[k]){
      return val[k];
    }
    if(j<col[k]){
      throw(invalid_argument("double& CSR::operator()(int, int) : bad indices\n"));
    }
  }
  throw(invalid_argument("double& CSR::operator()(int, int) : bad indices\n"));
}

CSR CSR::operator()(const vector<int>& I, const vector<int>& J) const{
  CSR A(static_cast<int>(I.size()), static_cast<int>(J.size()));
  int* col2 = new int[A.nr*A.nc];
  double* val2 = new double[A.nr*A.nc];
  int i, j;
  bool nullrow, nullprevrow;
  for(i=0;i<A.nr;i++){
    nullrow = true;
    for(j=0;j<A.nc;j++){
      if(this->operator()(I[i], J[j])!=0.0){
	nullrow = false;
	break;
      }
    }
    if(nullrow){
      if(i==0){
	A.row[i] = 0;
      }else{
	if(nullprevrow){
	  A.row[i] = A.row[i-1];
	}else{
	  A.row[i] = 1+A.nnz;
	}
      }
    }else{
      A.row[i] = A.nnz;
      col2[A.nnz] = j;
      val2[A.nnz] = this->operator()(I[i], J[j]);
      for(j++, A.nnz++;j<A.nc;j++){
	if(this->operator()(I[i], J[j])!=0.0){
	  col2[A.nnz] = j;
	  val2[(A.nnz)++] = this->operator()(I[i], J[j]);
	}
      }
    }
    nullprevrow = nullrow;
  }

  A.col = new int[A.nnz];
  A.val = new double[A.nnz];
  for(i=0;i<A.nnz;i++){
    A.col[i] = col2[i];
    A.val[i] = val2[i];
  }
  delete col2;
  delete val2;
  return A;
}

void CSR::MvProd(const vector<double>& x, vector<double>& b) const{
  assert((vector<double>::size_type)nc==x.size() && &x!=&b);

  int i, j;
  b.resize(0);
  b.reserve(nr);
  for(i=0;i+1<nr;i++){
    b.push_back(0.0);
    for(j=row[i];j<row[i+1];j++){
      b[i] += val[j]*x[col[j]];
    }
  }
  b.push_back(0.0);
  for(j=row[i];j<nnz;j++){//i==nr-1
    b[i] += val[j]*x[col[j]];
  }
}

void CSR::LUSolve(vector<double>& x, const vector<double>& b) const{
  assert(nc==nr && (vector<double>::size_type)nc==b.size());

  int i, j, k;
  DenseMatrix LU((DenseMatrix)(*this));

  //we compute the LU decomposition
  for(k=0;k<=nc-2;k++){
    for(i=1;i<nc-k;i++){
      LU(k+i, k) = LU(k+i, k)/LU(k, k);
      for(j=1;j<nc-k;j++){
	LU(k+i, k+j) = LU(k+i, k+j) - LU(k+i, k)*LU(k, j+k);
      }
    }
  }

  //we solve LX = b
  x = b;
  for(i = 1;i<nc;i++){
    for(j = 0;j<i;j++){
      x[i] -= LU(i, j)*x[j];
    }
  }

  //we solve Ux = X
  for(i = nc-1;i>=0;i--){
    for(j = i+1;j<nc;j++){
      x[i] -= LU(i, j)*x[j];
    }
    x[i] /= LU(i, i);
  }
}

CSR* CSR::clone() const{
  return new CSR{*this};
}

ostream& operator<<(ostream& os, const CSR& A){
  os << A.nr << " by " << A.nc << " Sparse Matrix (CSR):\nnumber of non-zeros: " << A.nnz;
  os << "\nrow:\n[ ";
  int i;
  if(A.nr){
    os << A.row[0];
  }
  for(i=1;i<A.nr;i++){
    os << ", " << A.row[i];
  }
  os << " ]\n(col, val):\n[ ";
  if(A.nnz){
    os << "( " << A.col[0] << ", " << A.val[0] << " )";
  }
  for(i=1;i<A.nnz;i++){
    os << ", ( " << A.col[i] << ", " << A.val[i] << " )";
  }
  os << " ]\n";
  return os;
}

CSR CSR::laplacian(int n){
  assert(n>=3);
  CSR A(n, n);
  //n>=3
  A.nnz = 3*n-2;
  A.col = new int[A.nnz];
  A.val = new double[A.nnz];
  int i;
  A.col[0] = 0;
  A.val[0] = 2.;
  A.col[1] = 1;
  A.val[1] = -1.;

  A.row[1] = 2;
  A.col[2] = 0;
  A.col[3] = 1;
  A.col[4] = 2;
  A.val[2] = -1.;
  A.val[3] = 2.;
  A.val[4] = -1.;
  for(i=2;i<n-1;i++){
    A.row[i] = A.row[i-1]+3;
    A.col[3*i-1] = i-1;
    A.col[3*i] = i;
    A.col[3*i+1] = i+1;
    A.val[3*i-1] = -1.;
    A.val[3*i] = 2.;
    A.val[3*i+1] = -1.;
  }

  A.row[n-1] = 3*(n-1)-1;
  A.col[3*(n-1)-1] = n-2;
  A.col[3*(n-1)] = n-1;
  A.val[3*(n-1)-1] = -1.;
  A.val[3*(n-1)] = 2.;

  return A;
}

CSR CSR::H(int n){
  assert(n>=4);
  CSR H(n, n);

  H.nnz = 3*n-4;

  H.col = new int[H.nnz];
  H.val = new double[H.nnz];

  H.col[0] = 0;
  H.val[0] = 1.0;
  H.col[1] = 2;
  H.val[1] = 1.0;

  H.row[1] = 2;
  H.col[2] = 1;
  H.val[2] = 2.0;
  H.col[3] = 3;
  H.val[3] = 1.0;
  for(int i=2;i<n-2;i++){
    H.row[i] = (i-2)*3 + 4;
    H.col[H.row[i]] = i-2;
    H.val[H.row[i]] = 1.0;
    H.col[H.row[i]+1] = i;
    H.val[H.row[i]+1] = static_cast<double>(i+1);
    H.col[H.row[i]+2] = i+2;
    H.val[H.row[i]+2] = 1.0;
  }

  //the last but one row
  H.row[n-2] = H.nnz-4;
  H.col[H.nnz-4] = n-4;
  H.val[H.nnz-4] = 1.0;
  H.col[H.nnz-3] = n-2;
  H.val[H.nnz-3] = static_cast<double>(n-1);

  //the last row
  H.row[n-1] = H.nnz-2;
  H.col[H.nnz-2] = n-3;
  H.val[H.nnz-2] = 1.0;
  H.col[H.nnz-1] = n-1;
  H.val[H.nnz-1] = static_cast<double>(n);

  return H;
}
