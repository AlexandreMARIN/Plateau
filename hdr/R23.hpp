#ifndef R23_HPP
#define R23_HPP

#include <cassert>


class N3 {
  int x, y, z;

public:
  N3() = default;
  N3(const N3&) = default;
  N3(N3&&) = default;
  N3(int xx, int yy, int zz) : x(xx), y(yy), z(zz){}

  N3& operator=(const N3&) = default;
  N3& operator=(N3&&) = default;

  int operator()(int) const;
  int& operator()(int);
};


class R2 {

  double x, y;

public:
  R2() = default;
  R2(const R2&) = default;
  R2(R2&&) = default;
  R2(double xx, double yy) : x{xx}, y{yy}{}

  R2& operator=(const R2&) = default;
  R2& operator=(R2&&) = default;
  R2 operator-() const;
  double operator()(int) const;
  double& operator()(int);
  friend R2 operator+(const R2&, const R2&);
  friend R2 operator-(const R2&, const R2&);
  friend R2 operator*(double, const R2&);
  friend double operator,(const R2&, const R2&);
  R2 ortho() const;
};


class R3 {

  double x, y, z;

public:
  R3() = default;
  R3(const R3&) = default;
  R3(R3&&) = default;
  R3(double xx, double yy, double zz) : x{xx}, y{yy}, z{zz}{}

  R3& operator=(const R3&) = default;
  R3& operator=(R3&&) = default;
  R3 operator-() const;
  double operator()(int) const;
  double& operator()(int);
  friend R3 operator+(const R3&, const R3&);
  friend R3 operator-(const R3&, const R3&);
  friend R3 operator*(double, const R3&);
  friend double operator,(const R3&, const R3&);
};


/*class N3*/
inline int N3::operator()(int i) const{
  assert(i>=0 && i<3);
  if(i==0){
    return x;
  }
  if(i==1){
    return y;
  }
  return z;
}

inline int& N3::operator()(int i){
  assert(i>=0 && i<3);
  if(i==0){
    return x;
  }
  if(i==1){
    return y;
  }
  return z;
}

/*class R2*/
inline R2 R2::operator-() const{
  return R2{-x, -y};
}

inline double R2::operator()(int i) const{
  assert(i>=0 && i<2);
  if(i==0){
    return x;
  }
  return y;
}

inline double& R2::operator()(int i){
  assert(i>=0 && i<2);
  if(i==0){
    return x;
  }
  return y;
}

inline R2 operator+(const R2& a, const R2& b){
  return R2{a.x + b.x, a.y + b.y};
}

inline R2 operator-(const R2& a, const R2& b){
  return R2{a.x - b.x, a.y - b.y};
}

inline R2 operator*(double a, const R2& u){
  return R2{a*u.x, a*u.y};
}

inline double operator,(const R2& a, const R2& b){
  return a.x*b.x + a.y*b.y;
}

inline R2 R2::ortho() const{
  return R2{-y, x};
}


/*class R3*/
inline R3 R3::operator-() const{
  return R3{-x, -y, -z};
}

inline double R3::operator()(int i) const{
  assert(i>=0 && i<3);
  if(i==0){
    return x;
  }
  if(i==1){
    return y;
  }
  return z;
}

inline double& R3::operator()(int i){
  assert(i>=0 && i<3);
  if(i==0){
    return x;
  }
  if(i==1){
    return y;
  }
  return z;
}

inline R3 operator+(const R3& a, const R3& b){
  return R3{a.x + b.x, a.y + b.y, a.z + b.z};
}

inline R3 operator-(const R3& a, const R3& b){
  return R3{a.x - b.x, a.y - b.y, a.z - b.z};
}

inline R3 operator*(double a, const R3& u){
  return R3{a*u.x, a*u.y, a*u.z};
}

inline double operator,(const R3& a, const R3& b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}


#endif
