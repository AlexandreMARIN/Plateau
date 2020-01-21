#include <limits>
#include <cassert>
#include <cmath>
#include <typeinfo>

#include "../include/iterative_solvers.hpp"

using namespace std;


/*Class IterSolver*/

double IterSolver::norm(const vector<double>& x, Norm norm_choice){
  double res = 0.0;

  switch(norm_choice){
  case Norm::EUCLIDIAN:
    for(double xi : x){
      res += xi*xi;
    }
    return sqrt(res);
  case Norm::INFTY:
    for(double xi : x){
      if(abs(xi)>res){
	res = abs(xi);
      }
    }
    return res;
  default://A-norm with A symmetric, positive
    if(typeid(*A_)==typeid(DenseMatrix)){
      const DenseMatrix* A2_ = dynamic_cast<const DenseMatrix*>(A_);
      int i, j;
      for(i = 0;(vector<double>::size_type)i<x.size();i++){
	for(j = 0;j<i;j++){
	  res += 2.0*(*A2_)(i, j)*x[i]*x[j];
	}
	res += (*A2_)(i, i)*x[i]*x[i];
      }
    }else{
      vector<double> Ax;
      A_->MvProd(x, Ax);
      for(int i=0;(vector<double>::size_type)i<x.size();i++){
	res += x[i]*Ax[i];
      }
    }
  }
  return sqrt(res);
}

double IterSolver::sq_norm(const vector<double>& x, Norm norm_choice){
  double res = 0.0;

  switch(norm_choice){
  case Norm::EUCLIDIAN:
    for(double xi : x){
      res += xi*xi;
    }
    return res;
  case Norm::INFTY:
    for(double xi : x){
      if(abs(xi)>res){
	res = abs(xi);
      }
    }
    return res*res;
  default://A-norm with A symmetric, positive
    if(typeid(*A_)==typeid(DenseMatrix)){
      const DenseMatrix* A2_ = dynamic_cast<const DenseMatrix*>(A_);
      int i, j;
      for(i = 0;(vector<double>::size_type)i<x.size();i++){
	for(j = 0;j<i;j++){
	  res += 2.0*(*A2_)(i, j)*x[i]*x[j];
	}
	res += (*A2_)(i, i)*x[i]*x[i];
      }
    }else{
      vector<double> Ax;
      A_->MvProd(x, Ax);
      for(int i=0;(vector<double>::size_type)i<x.size();i++){
	res += x[i]*Ax[i];
      }
    }

  }
  return res;
}

void IterSolver::set_A(const Matrix* A){
  A_ = A;
}

void IterSolver::set_b(const std::vector<double>* b){
  b_ = b;
}

void IterSolver::set_tol(double tol){
  if(tol>0.){
    tol_ = tol;
  }
}

void IterSolver::set_n_max(int n_max){
  if(n_max>0){
    n_max_ = n_max;
  }
}

void IterSolver::set_norm(Norm norm){
  norm_ = norm;
}

const std::vector<double>& IterSolver::get_x() const{
  return x_;
}

int IterSolver::get_niter() const{
  return niter_;
}

const vector<double>& IterSolver::get_resvec() const{
  return resvec_;
}

const std::vector<double>& IterSolver::get_r() const{
  return r_;
}

void IterSolver::solve(){
  assert(A_ && b_);
  this->init();
  double rnorm, bnorm=norm(*b_, norm_), limit=bnorm*tol_;

  resvec_.resize(0);
  resvec_.reserve(n_max_+1);
  niter_ = 0;
  x_.resize(b_->size());
  x_.assign(x_.size(), 0.0);
  r_ = *b_;

  while( (rnorm=norm(r_, norm_))>limit && niter_<n_max_ ){
    update_solution();
    update_resvec();
    resvec_.push_back((bnorm!=0.0)?rnorm/bnorm:numeric_limits<double>::quiet_NaN());
    niter_++;
  }

  resvec_.push_back((bnorm!=0.0)?rnorm/bnorm:numeric_limits<double>::quiet_NaN());

}


/*Class ConjGrad*/

ConjGrad ConjGrad::obj;

void ConjGrad::init(){
  p_ = *b_;//p_0 = r_0
  aux2.resize(p_.size());
}

void ConjGrad::update_solution(){
  aux = sq_norm(r_, Norm::EUCLIDIAN);
  alpha_ = aux/sq_norm(p_, Norm::NORM_A);
  for(int i=0;(vector<double>::size_type)i<x_.size();i++){
    x_[i] += alpha_*p_[i];
  }
}

void ConjGrad::update_resvec(){

  int i;
  A_->MvProd(p_, aux2);
  for(i=0;(vector<double>::size_type)i<r_.size();i++){
    r_[i] -= alpha_*aux2[i];
  }
  aux = sq_norm(r_, Norm::EUCLIDIAN)/aux;
  for(i=0;(vector<double>::size_type)i<r_.size();i++){
    p_[i] = r_[i] + aux*p_[i];
  }

}

ConjGrad& ConjGrad::getobj(){
  return obj;
}

ConjGrad& ConjGrad_ = ConjGrad::getobj();


/*Class MinRes*/

MinRes MinRes::obj;

void MinRes::init(){
  aux.resize(b_->size());
}

void MinRes::update_solution(){
  A_->MvProd(r_, aux);
  alpha_ = sq_norm(r_, Norm::NORM_A)/sq_norm(aux, Norm::EUCLIDIAN);
  for(int i=0;(vector<double>::size_type)i<x_.size();i++){
    x_[i] += alpha_*r_[i];
  }
}

void MinRes::update_resvec(){

  A_->MvProd(x_, aux);
  for(int i=0;(vector<double>::size_type)i<r_.size();i++){
    r_[i] = (*b_)[i] - aux[i];
  }

}

MinRes& MinRes::getobj(){
  return obj;
}

MinRes& MinRes_ = MinRes::getobj();
