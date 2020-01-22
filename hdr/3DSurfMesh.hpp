#ifndef 3DSURFMESH_HPP
#define 3DSURFMESH_HPP

#include <vector>
#include <string>

#include "R23.hpp"

class 3DSurfMesh {

public:
  3DSurfMesh() = default;
  3DSurfMesh(const 3DSurfMesh&) = default;
  3DSurfMesh(3DSurfMesh &&) = default;
  3DSurfMesh(std::string);


  3DSurfMesh& operator=(const 3DSurfMesh&) = default;
  3DSurfMesh& operator=(3DSurfMesh&&) = default;
  const std::vector<R3>& vert() const;
  const std::vector<bool>& isonboundary() const;
  const std::vector<N3>& tri() const;
  const std::vector<std::vector<int> >& vert_to_tri() const;

  void save(std::string) const;

private:
  std::vector<R3> vert;
  std::vector<bool> isonboundary;
  std::vector<N3> tri;
  std::vector<std::vector<int> > vert_to_tri;

};




#endif
