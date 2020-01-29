#ifndef PLATEAU_HPP
#define PLATEAU_HPP

#include <string>
#include "SurfMesh3D.hpp"

class Plateau {

public:
  Plateau() = default;
  Plateau(const Plateau&) = default;
  Plateau(Plateau&&) = default;
  Plateau(const std::string&, double, int);

  Plateau& operator=(const Plateau&) = default;
  Plateau& operator=(Plateau&&) = default;

  void solve();
  void exportGnuplot(const std::string&);
  int get_iter();

private:

  SurfMesh3D surf;
  double tol_;
  int maxnbiter_;
  int iter_;

};

inline int Plateau::get_iter(){
  return iter_;
}

#endif
