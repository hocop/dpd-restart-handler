//
// Created by rusl on 3/21/17.
//

#ifndef LPD_SYSTEM_H
#define LPD_SYSTEM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include "vect3/vect3.h"
#include "Conf.h"

#define TABSIZE 8

using namespace std;

struct Restart {
	Restart();
	Restart(Conf const &conf);
	Restart(Restart const &r);
	
	int natoms; // number of atoms
	int nbonds; // number of bonds
	int nangles; // number of angles
	int nchains; // number of chains (sequence where each atom has 2 bonds)
	int nnodes; // number of nodes (when 1 atom has >2 bonds)
	
	flt density;
	
	vector<vect3> x;
	vector<int> valency, type, flag;
	vector<int> a1, a2; // list of atoms that are bonded
	vect3 box;
	
	vector<vector<int> > chain; // sequences of particles joined in one chain
	vector<int> center; // centers of nodes
	vector<vector<int> > shell; // list of all shell particles
	vector<vector<int> > bond; // list of bonds inside node which we want to simulate (others are computed somewhere else)
	
	void distPBC(vect3 &x) const;
	void absPBC_while(vect3 &x) const;
	void absPBC_if(vect3 &x) const;
	
	void add(Restart const& r, vect3 r0, flt a, flt b, flt c, flt scale=1.f); //r0 - pos; a,b,c - Euler angles
};

Restart importDpdNano(string fname, int withflag, bool withTopology = true); // dpd nano dat file
Restart importMol(string fname, Conf conf); // mol file
Restart importSmall(string fname); // smallest possible format: ~1/10 of gav. dat file size
Restart importNormal(string fname); // format that can be read with eyes: ~1/3 of gav. dat file size
Restart importLmp(string fname);

void exportGav(Restart &rst, string fname, int withflag);
void exportSmall(Restart &rst, string fname);
void exportNormal(Restart &rst, string fname);

void topology(Restart &r); // make chains out of bonds
void topologyReverse(Restart &r); // make bonds out of chains

#endif //LPD_SYSTEM_H
