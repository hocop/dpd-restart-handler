//
// Created by rusl on 4/11/17.
//

#ifndef ANALYSE_ANALYSER_H
#define ANALYSE_ANALYSER_H

#include "Restart.h"

// Systems
Restart selectByType(Restart const &r, vector<int> &type); // get system with only particles of certain types
// you can add one system to another using method Restart::add

// Common math operations
void meanAndStdev(vector<flt> value, flt &mean, flt &stdev); // get mean and stdev of this array

// Chains
vect3 centerOfMassOfChain(Restart const &r, int i); // i is index of chain
void countEEDist(Restart const &r, flt &mean, flt &stdev); // end-to-end distance; mean and stdev

// Liquid crystals
vector<flt> countNematicOrderBonds(Restart const &r, int rmin, int rmax);
vector<flt> countNematicOrderChains(Restart const &r, int rmin, int rmax);
vect3 countDirectorChains(Restart const &r);
vect3 countDirectorBonds(Restart const &r); // not implemented
flt countSmecticOrder(Restart const &r, flt periodRec);


#endif //ANALYSE_ANALYSER_H
