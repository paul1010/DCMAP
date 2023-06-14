#include <iostream>
#include <armadillo>
#include <algorithm>
#include <vector>

using namespace std;
using namespace arma;

#include <dcmap.h>

int main(int argc, char** argv)
  {
   
    // ---------------------------- Simple Example  ---------------------------- //
    std::vector<std::string> names;
    names = { "A", "B", "C", "D", "E", "F", "G" };
    uword numnodes = names.size();

    mat cptparchild(numnodes, numnodes, fill::value(datum::nan));
    cptparchild(4, 2) = 0; // E parent is C
    cptparchild(3, span(0, 1)) = mat(1,2,fill::value(0)); // D parent is A, B
    cptparchild(5, 0) = 0; // F parent is 0
    cptparchild(6, span(3, 4)) = mat(1,2,fill::value(0)); // G parent is D,E
        
    uvec numstates(numnodes, fill::value(2));
    uvec leaves = { 5,6 };
    uvec roots = { 0,1,2 };
    
    // Print simple example
    //for (int i = 0; i < names.size(); ++i)
    //    cout << names[i] << " ";
    //cout << "numnodes = " << numnodes << endl;
    //cptparchild.print("cptparchild = ");
    //numstates.print("numstates = ");
    //leaves.print("leaves = ");
    //roots.print("roots = ");

    // Constants
    static int maxiter = 1000;

    // Run
    dcmap D(numnodes, numstates, leaves, roots, cptparchild,maxiter);
    D.search(1);
    
    D.debugdat.save("C:/Work/debugdat.csv", csv_ascii);

  return 0;
  }





  
