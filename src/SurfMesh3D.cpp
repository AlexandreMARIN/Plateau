#include "../hdr/SurfMesh3D.hpp"

#include <fstream>
#include <set>
#include <iostream>

using namespace std;


SurfMesh3D::SurfMesh3D(const string& filename) : vert_(0), isonboundary_(0), tri_(0), vert_to_tri_(0){

  ifstream file(filename);

  string buf;
  int vertnb, trinb;
  const int maxnbneighbor = 7;
  string tag1{"Vertices"}, tag2{"Triangles"};

  do{
    getline(file, buf);
  }while(buf!=tag1);

  file >> vertnb;
  vert_.reserve(vertnb);

  double x, y, z;
  int label;
  for(int i=0;i<vertnb;i++){
    file >> x >> y >> z >> label;
    vert_.push_back(R3{x, y, z});
  }

  do{
    getline(file, buf);
  }while(buf!=tag2);

  file >> trinb;
  tri_.reserve(trinb);
  vert_to_tri_.resize(vertnb);

  int ind[3];
  set<pair<int, int> > edges;
  pair<int, int> edge;
  for(int k=0;k<trinb;k++){
    file >> ind[0] >> ind[1] >> ind[2] >> label;
    ind[0]--;
    ind[1]--;
    ind[2]--;
    tri_.push_back(N3{ind[0], ind[1], ind[2]});
    for(int l=0;l<3;l++){
      if(vert_to_tri_[ind[l]].empty()){
	vert_to_tri_[ind[l]].reserve(maxnbneighbor);
      }
      vert_to_tri_[ind[l]].push_back(k);
      if(ind[l] < ind[(l+1)%3]){
	edge = make_pair(ind[l], ind[(l+1)%3]);
      }else{
	edge = make_pair(ind[(l+1)%3], ind[l]);
      }
      if(edges.find(edge) == edges.end()){
	edges.insert(edge);
      }else{
	edges.erase(edge);
      }
    }
  }

  isonboundary_.resize(vertnb);
  for(pair<int, int> e : edges){
    isonboundary_[e.first] = true;
    isonboundary_[e.second] = true;
  }


}

void SurfMesh3D::set_z_coord(const vector<int>& index, const vector<double>& z_coord){
  for(int i=0;static_cast<vector<int>::size_type>(i)<index.size();i++){
    vert_[index[i]](2) = z_coord[i];
  }
}

void SurfMesh3D::exportGnuplot(const string& filename) const{

  ofstream file{filename};

  for(N3 n : tri_){
    for(int l=0;l<3;l++){
      file << vert_[n(l)](0) << " " << vert_[n(l)](1) << " " << vert_[n(l)](2) << "\n";
    }
    file << vert_[n(0)](0) << " " << vert_[n(0)](1) << " " << vert_[n(0)](2) << "\n\n\n";
  }

}

void SurfMesh3D::save(const string& filename) const{

  //in .mesh files, indices of vertices are positive !
  ofstream file{filename};

  file << "MeshVersionFormatted 2\n\nDimension 3\n\n";
  file << "\nVertices\n" << vert_.size() << "\n";

  for(const R3& v : vert_){
    file << v(0) << " " << v(1) << " "<< v(2) << " 0\n";
  }

  file << "\nTriangles\n" << tri_.size() << "\n";

  set<pair<int, int> > edges;
  pair<int, int> edge;
  for(const N3& n : tri_){
    file << n(0)+1 << " " << n(1)+1 << " " << n(2)+1 << " 0\n";//label: 0
    for(int l=0;l<3;l++){
      if(n(l%3) < n((l+1)%3)){
	edge = make_pair(n(l)+1, n((l+1)%3)+1);
      }else{
	edge = make_pair(n((l+1)%3)+1, n(l)+1);
      }
      if(edges.find(edge) == edges.end()){
	edges.insert(edge);
      }else{
	edges.erase(edge);
      }
    }
  }

  file << "\nEdges\n" << edges.size() << "\n";
  for(const auto& e : edges){
    file << e.first << " " << e.second << " 0\n";//label: 0
  }

  file << "\nEnd";
}
