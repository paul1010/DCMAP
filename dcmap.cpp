#include <iostream>
#include <armadillo>
#include <algorithm>
#include <vector>

using namespace std;
using namespace arma;

#include <dcmap.h>

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

// NOTE: the C++11 "auto" keyword is not recommended for use with Armadillo objects and functions

// NOTE: swap branch i and elements j as sparse matrices work faster (elements in column contiguous in memory)

int main(int argc, char** argv)
  {
    //cout << "Armadillo version: " << arma_version::as_string() << endl;
    // ---------------------------- Simple Example  ---------------------------- //
    std::vector<std::string> names;
    names = { "A", "B", "C", "D", "E", "F", "G" };
    uword numnodes = names.size();
    uvec _numnodes;
    _numnodes.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/numnodes.csv", csv_ascii);
    //_numnodes.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/numnodes2.csv", csv_ascii);
    numnodes = _numnodes(0);

    mat cptparchild(numnodes, numnodes, fill::value(datum::nan));
    //cptparchild(4, 2) = 0; // E parent is C
    //cptparchild(3, span(0, 1)) = mat(1,2,fill::value(0)); // D parent is A, B
    //cptparchild(5, 0) = 0; // F parent is 0
    //cptparchild(6, span(3, 4)) = mat(1,2,fill::value(0)); // G parent is D,E
        
    cptparchild.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/cptparchild.csv", csv_ascii);
    //cptparchild.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/cptparchild2.csv", csv_ascii);
    uvec numstates(numnodes, fill::value(2));
    numstates.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/numstates.csv", csv_ascii);
    //numstates.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/numstates2.csv", csv_ascii);
    uvec leaves = { 5,6 };
    leaves.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/leaves.csv", csv_ascii);
    //leaves.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/leaves2.csv", csv_ascii);
    uvec roots = { 0,1,2 };
    roots.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/roots.csv", csv_ascii);
    //roots.load("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/roots2.csv", csv_ascii);
    
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
    
    D.debugdat.save("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/debugdat.csv", csv_ascii);
    //(D.Qdash.exportq("---noprint")).save("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/Qdash.csv", csv_ascii);
    //(D.Q.exportq("---noprint")).save("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/Q.csv", csv_ascii);

    //uvec qqinds(D.Q.qinds);
    //qqinds.save("C:/Work/Seagrass Model 2/DBPCN_GPU/Cpp_Test/qinds.csv", csv_ascii);
    


  return 0;
  }





  
