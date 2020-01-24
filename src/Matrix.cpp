#include <chrono>
#include <random>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iterator>

#include "../hdr/Matrix.hpp"
#include "../hdr/iterative_solvers.hpp"

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

void DenseMatrix::resize(int r, int c){
  if(r*c != nr*nc){
    delete[] val;
    val = new double[r*c];
  }
  nr = r;
  nc = c;
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

CSR::CSR(int r, int c) : AtomicMatrix(r, c), nnz(0), row(0), colval{}{
  row.resize(r, colval.end());
}

CSR::CSR(const CSR& A) : AtomicMatrix(A.nr, A.nc), nnz(A.nnz), row(0), colval(A.colval){
  row.resize(nr, colval.end());
  auto itA = A.colval.begin();
  auto it = colval.begin();

  for(int i=0;i<nr;i++){
    if(A.row[i]==A.colval.end()){
      row[i] = colval.end();
      continue;
    }
    while(itA!=A.row[i]){
      it++;
      itA++;
    }
    row[i] = it;
  }
}

CSR& CSR::operator=(const CSR& A){
  nr = A.nr;
  nc = A.nc;
  nnz = A.nnz;
  colval = A.colval;
  row.resize(nr, colval.end());
  auto itA = A.colval.begin();
  auto it = colval.begin();

  for(int i=0;i<nr;i++){
    if(A.row[i]==A.colval.end()){
      row[i] = colval.end();
      continue;
    }
    while(itA!=A.row[i]){
      it++;
      itA++;
    }
    row[i] = it;
  }
  return *this;
}

double CSR::operator()(int i, int j) const{
  assert( (i>=0) && (i<nr) && (j>=0) && (j<nc) );

  if(row[i]==colval.end()){
    return 0.0;
  }
  int i2 = i+1;
  while( (i2<nr) && (row[i2]==colval.end()) ){
    i2++;
  }
  if(i2==nr){
    for(auto it=row[i];it!=colval.end();it++){
      if(j == it->first){
	return it->second;
      }
    }
    return 0.0;
  }

  for(auto it=row[i];it!=row[i2];it++){
    if(j == it->first){
      return it->second;
    }
    if(j < it->first){
      break;
    }
  }
  return 0.0;
}

double& CSR::operator()(int i, int j){
  assert( (i>=0) && (i<nr) && (j>=0) && (j<nc) );

  int i2;
  if(row[i]==colval.end()){
    i2 = i+1;
    while( (i2<nr) && (row[i2]==colval.end()) ){
      i2++;
    }
    if(i2<nr){
      nnz++;
      row[i] = colval.emplace(row[i2], j, 0.0);
      return row[i]->second;
    }
    //we insert the last line
    nnz++;
    row[i] = colval.emplace(colval.end(), j, 0.0);
    return row[i]->second;
  }

  i2 = i+1;
  //the i-th row exists
  while( (i2<nr) && (row[i2]==colval.end()) ){
    i2++;
  }
  if(i2==nr){//as a matter of fact, the i-th row is the last non zero line
    for(auto it=row[i];it!=colval.end();it++){
      if(j == it->first){
	return it->second;
      }
      if(j < it->first){
	nnz++;
	colval.emplace(it, j, 0.0);
	if(row[i]==it){
	  row[i]--;
	}
	it--;
	return it->second;
      }
    }
    //we add the last value
    nnz++;
    colval.emplace_back(j, 0.0);
    auto it = colval.end();
    it--;
    return it->second;
  }

  for(auto it=row[i];it!=row[i2];it++){
    if(j == it->first){
      return it->second;
    }
    if(j < it->first){
      nnz++;
      colval.emplace(it, j, 0.0);
      if(row[i]==it){
	row[i]--;
      }
      it--;
      return it->second;
    }
  }

  //we insert the last non zero element of the i-th row
  auto it = colval.emplace(row[i2], j, 0.0);
  nnz++;
  return it->second;
}


void CSR::MvProd(const vector<double>& x, vector<double>& b) const{
  assert(static_cast<vector<double>::size_type>(nc)==x.size() && &x!=&b);
  b.reserve(nr);
  b.resize(0);
  for(int i=0;i<nr-1;i++){
    b.push_back(0.0);
    for(auto it=row[i];it!=row[i+1];it++){
      b[i] += (it->second)*x[it->first];
    }
  }
  b.push_back(0.0);
  for(auto it=row[nr-1];it!=colval.end();it++){
    b[nr-1] += (it->second)*x[it->first];
  }
}


CSR* CSR::clone() const{
  return new CSR{*this};
}

ostream& operator<<(ostream& os, const CSR& A){
  os << A.nr << " by " << A.nc << " CSR matrix:\n";
  os << "(dereferenced) row:\n[ ";
  if(A.row.front()==A.colval.end()){
    os << "[]";
  }else{
    os << "[" << A.row.front()->first << ", " << A.row.front()->second << "]";
  }
  auto it = A.row.begin();
  for(it++;it!=A.row.end();it++){
    if(*it==A.colval.end()){
      os << ", []";
      continue;
    }
    os << ", [" << (*it)->first << ", " << (*it)->second << "]";
  }
  os << " ]\n";
  os << "colval:\n[ [" << A.colval.front().first << ", " << A.colval.front().second << "]";
  auto itcv = A.colval.begin();
  itcv++;
  for(;itcv!=A.colval.end();itcv++){
    os << ", [" << itcv->first << ", " << itcv->second << "]";
  }
  os << " ]\n";
  return os;
}


