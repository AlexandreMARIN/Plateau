#ifndef PLATEAU_HPP
#define PLATEAU_HPP

#include <string>
#include "SurfMesh3D.hpp"

/*
  That class enables us to solve Plateau's problem.
*/
class Plateau {

public:
  Plateau() = default;
  Plateau(const Plateau&) = default;
  Plateau(Plateau&&) = default;

  /*
    To that constructor, we must give:
    - a .mesh file
    - a tolerance tol_
    - a maximum number of iterations
    This class computes a sequence of n-tuples x_i and the algorithm ends when the infinity norm of (x_(i+1) - x(i)) is less than tol.
  */
  Plateau(const std::string&, double, int);

  Plateau& operator=(const Plateau&) = default;
  Plateau& operator=(Plateau&&) = default;

  //This method solves the problem and modifies the z-coordinates of each internal point of the mesh
  void solve();

  //That method writes a .mesh file to save the current mesh
  void exportGnuplot(const std::string&);

  int get_iter();//gives the number of iterations used during the previous call to solve()

private:

  SurfMesh3D surf;//the mesh we transform
  double tol_;//tolerance
  int maxnbiter_;//maximum number of iterations
  int iter_;//number of used iterations

};

inline int Plateau::get_iter(){
  return iter_;
}

#endif
