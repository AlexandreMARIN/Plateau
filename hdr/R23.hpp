#ifndef R23_HPP
#define R23_HPP

class R2 {

  double x, y;

public:
  R2() = default;
  R2(const R2&) = default;
  R2(R2&&) = default;
  R2(double xx, double yy) : x{xx}, y{yy}{}

  const R2& operator=(const R2&) = default;
  const R2& operator=(R2&&) = default;
  const R2& operator-() const;
  double operator()(int) const;
  double& operator()(int);
  friend R2 operator+(const R2&, const R2&);
  friend R2 operator-(const R2&, const R2&);
  double operator,(const R2&, const R2&);
  R2 ortho() const;
};


class R3 {

  double x, y, z;

public:
  R3() = default;
  R3(const R3&) = default;
  R3(R3&&) = default;
  R3(double xx, double yy, double zz) : x{xx}, y{yy}, z{zz}{}

  const R3& operator=(const R3&) = default;
  const R3& operator=(R3&&) = default;
  const R3& operator-() const;
  double operator()(int) const;
  double& operator()(int);
  friend R3 operator+(const R3&, const R3&);
  friend R3 operator-(const R3&, const R3&);
  double operator,(const R3&, const R3&);
};



inline const R2& operator-() const{
  x = -x;
  y = -y;
  return *this;
}

inline double operator()(int i) const{
  cassert(i>=0 && i<2);
  if(i==0){
    return x;
  }
  return y;
}

inline double& operator()(int i){
  cassert(i>=0 && i<2);
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

inline double operator,(const R2& a, const R2& b){
  return a.x*b.x + a.y*b.y;
}

inline R2 ortho() const{
  return R2{-y, x};
}



inline const R3& operator-() const{
  x = -x;
  y = -y;
  z = -z;
  return *this;
}

inline double operator()(int i) const{
  cassert(i>=0 && i<3);
  if(i==0){
    return x;
  }
  if(i==1){
    return y;
  }
  return z;
}

inline double& operator()(int i){
  cassert(i>=0 && i<3);
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

inline double operator,(const R3& a, const R3& b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}


#endif
