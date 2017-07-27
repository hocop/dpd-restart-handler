#include "Analyser.h"

// compile command:
// g++ EXAMPLE.cpp Analyser.cpp Restart.cpp Conf.cpp vect3/vect3.cpp -o example

// This is an example of how to use this lib
// It loads file restart.dat and computes number of chains in it
// Then it saves file as restart.sdat in another format

int main() {
	Restart restart = importGav("restart.dat"); // load file
	cout << "n chains: " << restart.chain.size() << endl; // print number of chains
	exportSmall(restart, "restart.sdat"); // save in new format
	return 0;
}
