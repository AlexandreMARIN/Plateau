#ifndef ITERATIVE_SOLVERS
#define ITERATIVE_SOLVERS

#include "Matrix.hpp"

/*The class IterSolver represents iterative solvers.*/
/*
  That class enables us to solve the linear system A_*x_ = b_.
*/

class IterSolver{
public:
  enum class Norm {INFTY, EUCLIDIAN, NORM_A};//for the choice of the norm
protected:
  //input
  const Matrix* A_;
  const std::vector<double>* b_;
  double tol_;//tolerance
  int n_max_;//maximal number of iterations

  //output
  std::vector<double> x_;//solution
  int niter_;//number of iterations
  std::vector<double> resvec_;//norms of the residual vector for each iterations

  std::vector<double> r_;//residual vector

  Norm norm_;//norm used for computing the quantity |r|/|b|

  IterSolver() = default;
  IterSolver(const IterSolver&) = delete;
  IterSolver(IterSolver&&) = delete;
  virtual ~IterSolver() = default;
  IterSolver& operator=(const IterSolver&) = delete;
  IterSolver& operator=(IterSolver&&) = delete;

  double norm(const std::vector<double>&, Norm);//returns the norm of the first parameter
  double sq_norm(const std::vector<double>&, Norm);//the squared norm of the fist parameter

  virtual void init() = 0;//function called at the beginning of the method solve()
  virtual void update_solution() = 0;//computes the new x_
  virtual void update_resvec() = 0;//computes the new r_

public:

  //setters
  void set_A(const Matrix*);
  void set_b(const std::vector<double>*);
  void set_tol(double);
  void set_n_max(int);
  void set_norm(Norm norm);

  //getters
  const std::vector<double>& get_x() const;
  int get_niter() const;
  const std::vector<double>& get_resvec() const;
  const std::vector<double>& get_r() const;

  void solve();

};


/***************Conjugate Gradient*****************/
class ConjGrad final: public IterSolver{

  //for computing
  double alpha_;
  double aux;
  std::vector<double> p_, aux2;

  static ConjGrad obj;

  ConjGrad() = default;
  ConjGrad(const ConjGrad&) = delete;
  ConjGrad(ConjGrad&&) = delete;
  ~ConjGrad() = default;
  ConjGrad& operator=(const ConjGrad&) = delete;
  ConjGrad& operator=(ConjGrad&&) = delete;

  void init() override;
  void update_solution() override;
  void update_resvec() override;

public:

  static ConjGrad& getobj();

};

extern ConjGrad& ConjGrad_;


/***********Method MinRes***************/
class MinRes final: public IterSolver{

  //for computing
  double alpha_;
  std::vector<double> aux;

  static MinRes obj;

  MinRes() = default;
  MinRes(const MinRes&) = delete;
  MinRes(MinRes&&) = delete;
  ~MinRes() = default;
  MinRes& operator=(const MinRes&) = delete;
  MinRes& operator=(MinRes&&) = delete;

  void init() override;
  void update_solution() override;
  void update_resvec() override;

public:

  static MinRes& getobj();

};

extern MinRes& MinRes_;


#endif
