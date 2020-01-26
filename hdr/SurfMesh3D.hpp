#ifndef SURFMESH3D_HPP
#define SURFMESH3D_HPP

#include <vector>
#include <string>

#include "R23.hpp"

class SurfMesh3D {

public:
  SurfMesh3D() = default;
  SurfMesh3D(const SurfMesh3D&) = default;
  SurfMesh3D(SurfMesh3D &&) = default;
  SurfMesh3D(const std::string&);


  SurfMesh3D& operator=(const SurfMesh3D&) = default;
  SurfMesh3D& operator=(SurfMesh3D&&) = default;
  //void set_z_coord_vert(const std::vector<double>&);
  const std::vector<R3>& vert() const;
  const std::vector<bool>& isonboundary() const;
  const std::vector<N3>& tri() const;
  const std::vector<std::vector<int> >& vert_to_tri() const;

  //void save(std::string) const;

private:
  std::vector<R3> vert_;
  std::vector<bool> isonboundary_;
  std::vector<N3> tri_;
  std::vector<std::vector<int> > vert_to_tri_;

};


inline const std::vector<R3>& SurfMesh3D::vert() const{
  return vert_;
}

inline const std::vector<bool>& SurfMesh3D::isonboundary() const{
  return isonboundary_;
}

inline const std::vector<N3>& SurfMesh3D::tri() const{
  return tri_;
}

inline const std::vector<std::vector<int> >& SurfMesh3D::vert_to_tri() const{
  return vert_to_tri_;
}


#endif
