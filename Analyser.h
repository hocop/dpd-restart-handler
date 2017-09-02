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

// Chains low-level operations
vect3 centerOfMassOfChain(Restart const &r, int i); // i is index of chain
void countEEDist(Restart const &r, flt &mean, flt &stdev); // end-to-end distance; mean and stdev

// Common termodynamics
enum Units {ATOMS, CHAINS, NODES};
enum Space {LINE, PLANE, VOLUME};
vector<vect3> countDisplacement(Restart &r1, Restart &r2, Units units, Space space=VOLUME, vect3 dir=vect3());

// Liquid crystals
vector<flt> countNematicOrderBonds(Restart const &r, int rmin, int rmax);
vector<flt> countNematicOrderChains(Restart const &r, int rmin, int rmax);
flt countNematicOrderChainsMaxR(Restart const &r);
vect3 countDirectorChains(Restart const &r);
vect3 countDirectorBonds(Restart const &r); // not implemented
flt countSmecticOrder(Restart const &r, flt periodRec);

//Polymer
void countPDI (Restart const &r ,flt &avg, flt &PDI) ;

#endif //ANALYSE_ANALYSER_H
