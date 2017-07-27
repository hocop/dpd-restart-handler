#include "Analyser.h"

// compile command:
// g++ EXAMPLE.cpp Analyser.cpp Restart.cpp Conf.cpp vect3/vect3.cpp -o example

// This is to demonstrate usage of the lib

int main() {
	// load file
	Restart restart = importGav("restart.dat");
	
	// print number of chains
	cout << "n chains: " << restart.chain.size() << endl;
	
	// compute order parameter (if it is nematic then it will be 1, else 0)
	float order = countNematicOrderChainsMaxR(restart);
	cout << "order par: " << order << endl;
	
	// save new format
	exportSmall(restart, "restart.sdat");
	return 0;
}
