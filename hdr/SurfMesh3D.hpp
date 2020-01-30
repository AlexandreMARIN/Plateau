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

  /*
    This constructor reads a file with extension .mesh
  */
  SurfMesh3D(const std::string&);


  SurfMesh3D& operator=(const SurfMesh3D&) = default;
  SurfMesh3D& operator=(SurfMesh3D&&) = default;

  /* That function sets z-coordinates of nodes whose indices are given in the first argument */
  void set_z_coord(const std::vector<int>&, const std::vector<double>&);

  /* Some methods to access data members */
  const std::vector<R3>& vert() const;
  const std::vector<bool>& isonboundary() const;
  const std::vector<N3>& tri() const;
  const std::vector<std::vector<int> >& vert_to_tri() const;

  /* That method writes a file to represent the mesh.
     Then the created file can be used by gnuplot.
 */
  void exportGnuplot(const std::string&) const;

private:
  std::vector<R3> vert_;//vertices
  std::vector<bool> isonboundary_;//isonboundary[i] : is the i-th node in vert_ on the boundary ?
  std::vector<N3> tri_;//triangles
  std::vector<std::vector<int> > vert_to_tri_;//the i-th cell gives all the indices in tri_ of the triangles which contain the i-th node in vert_

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
