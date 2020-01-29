#include "../hdr/Plateau.hpp"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]){


  //we extract some arguments
  if(argc != 5){
    cerr << "There are not enough parameters: four arguments are expected.\nUSAGE: ./solvePlateau file_in tolerance maxnbiter file_out\n";
    return 1;
  }

  string infilename{argv[1]}, arg2{argv[2]}, arg3{argv[3]}, outfilename{argv[4]};

  int maxnbiter{stoi(arg3)};
  double tol{stod(arg2)};

  //we solve Plateau's problem and we export results
  Plateau problem{infilename, tol, maxnbiter};
  problem.solve();
  problem.exportGnuplot(outfilename);

  cout << "\nResults have been acquired in " << problem.get_iter() << " iterations.\n\n";

  return 0;
}
