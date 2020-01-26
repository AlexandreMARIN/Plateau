#include "../hdr/Plateau.hpp"
#include "../hdr/iterative_solvers.hpp"

#include <cmath>
#include <set>

using namespace std;

Plateau::Plateau(const string& filename, double eps, int maxnbiter) : surf{filename}, tol_{eps}, maxnbiter_{maxnbiter}, iter_{0}{}


void Plateau::solve(){

  //we retrieve information from the mesh
  const vector<R3>& vert = surf.vert();
  const vector<bool>& isonboundary = surf.isonboundary();
  const vector<N3>& tri = surf.tri();
  const vector<vector<int> >& vert_to_tri = surf.vert_to_tri();

  vector<double> w_(vert.size()), u_;//it will contain the solutions

  //we prepare two arrays to represent a 1-1 map which enables us to associate an index of an internal node to an index in the range [0, number of internal nodes - 1]
  vector<int> v2i, i2v;

  v2i.reserve(vert.size());
  i2v.reserve(vert.size());
  u_.reserve(vert.size());

  /* we fill u_ with the z-coordinates of all the nodes on the boundary in the initial mesh */
  for(int i=0; static_cast<vector<bool>::size_type>(i)<isonboundary.size() ;i++){
    if(isonboundary[i]){
      v2i.push_back(i2v.size());
      i2v.push_back(i);
      u_.push_back(vert[i](2));
    }else{
      v2i.push_back(-1);
    }
  }

  //we allocate enough space for the vector b_ and we declare some variables
  CSR A_(static_cast<int>(i2v.size()), static_cast<int>(i2v.size()));//a matrix in the CSR format
  vector<double> b_(i2v.size());

  int li;//see below
  double diff, diff2;//for the stop criterion
  vector<double> normN_K(tri.size());//for each triangle K, it contains the norm of the normal vector when w_ gives the z-coordinates
  double N_K3;//the z-coordinate of the normal vector
  R2 N_K12, P01, P02, DN_K12i, DN_K12j;//we use the same notations as in the wording of the problem
  set<int> colind;//for each row index i, it will contain the indices j such that A_(i, j) is not zero

  iter_ = 0;

  do{

    iter_++;
    w_ = u_;

    //we build the linear system

    //for each triangle, we compute the norm of the normal vector
    for(int k=0;static_cast<vector<int>::size_type>(k)<tri.size();k++){

      N_K12 = R2{0., 0.};
      for(int l=0;l<3;l++){
	if(isonboundary[tri[k](l)]){
	  N_K12 = N_K12 + vert[tri[k](l)](2) * R2{vert[tri[k]((l+1)%3)](0) - vert[tri[k]((l+2)%3)](0), vert[tri[k]((l+1)%3)](1) - vert[tri[k]((l+2)%3)](1)};
	}else{
	  N_K12 = N_K12 + w_[v2i[tri[k](l)]] * R2{vert[tri[k]((l+1)%3)](0) - vert[tri[k]((l+2)%3)](0), vert[tri[k]((l+1)%3)](1) - vert[tri[k]((l+2)%3)](1)};
	}
      }

      P01 = R2{vert[tri[k](1)](1) - vert[tri[k](0)](1), vert[tri[k](1)](2) - vert[tri[k](0)](2)};
      P02 = R2{vert[tri[k](2)](1) - vert[tri[k](0)](1), vert[tri[k](2)](2) - vert[tri[k](0)](2)};
      N_K3 = P01(0)*P02(1) - P01(1)*P02(0);

      normN_K[k] = sqrt( N_K12(0)*N_K12(0) + N_K12(1)*N_K12(1) + N_K3*N_K3 );

    }


    //we assembly the matrix A_ and the vector b_
    for(int i = 0 ; i < A_.get_nr() ; i++){

      A_(i, i) = 0.0;
      b_[i] = 0.0;
      colind.clear();
      for(int tri_ind : vert_to_tri[i2v[i]]){

	//we find the i-th internal node and we compute a part of A_(i, i)
	for(int l=0;l<3;l++){
	  if(tri[tri_ind](l) == i2v[i]){
	    li = l;
	    DN_K12i = R2{vert[tri[tri_ind]((l+1)%3)](0) - vert[tri[tri_ind]((l+2)%3)](0), vert[tri[tri_ind]((l+1)%3)](1) - vert[tri[tri_ind]((l+2)%3)](1)};
	    A_(i, i) += (DN_K12i, DN_K12i) / normN_K[tri_ind];
	    break;
	  }
	}

	//we compute a part of A_(i, j) and b_[i]
	for(int l=0;l<3;l++){

	  if(l==li){
	    continue;
	  }
	  if(isonboundary[tri[tri_ind](l)]){
	    N_K12 = vert[tri[tri_ind](l)](2) * R2{vert[tri[tri_ind]((l+1)%3)](0) - vert[tri[tri_ind]((l+2)%3)](0), vert[tri[tri_ind]((l+1)%3)](1) - vert[tri[tri_ind]((l+2)%3)](1)};
	    b_[i] += (DN_K12i, N_K12) / normN_K[tri_ind];
	  }else{
	    if(v2i[tri[tri_ind](l)] < i){
	      if( colind.find(v2i[tri[tri_ind](l)]) == colind.end() ){
		A_(i, v2i[tri[tri_ind](l)]) = 0.;
		colind.insert(v2i[tri[tri_ind](l)]);
	      }
	      DN_K12j = R2{vert[tri[tri_ind]((l+1)%3)](0) - vert[tri[tri_ind]((l+2)%3)](0), vert[tri[tri_ind]((l+1)%3)](1) - vert[tri[tri_ind]((l+2)%3)](1)};
	      A_(i, v2i[tri[tri_ind](l)]) += (DN_K12i, DN_K12j) / normN_K[tri_ind];
	    }
	  }

	}

      }

      b_[i] = -b_[i];
      //A_ is symmetric
      for(int j : colind){
	A_(j, i) = A_(i, j);
      }

    }

    //we solve the linear system
    ConjGrad_.set_n_max(A_.get_nr());
    A_.cg(u_, b_);

    //we compute a norm of (u_ - w_) to check a condition of the stop criterion
    diff = abs(u_[0] - w_[0]);
    for(int i=1;i<A_.get_nr();i++){
      diff2 = abs(u_[i] - w_[i]);
      if(diff < diff2){
	diff = diff2;
      }
    }

  }while( (iter_ < maxnbiter_) && (diff > tol_) );

}
