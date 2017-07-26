//
// Created by rusl on 3/21/17.
//

#ifndef LPD_CONF_H
#define LPD_CONF_H

#include <vector>
#include "vect3/vect3.h"

enum RUN_TYPE {re, co};

struct Conf {
	RUN_TYPE run_type;
	int species_number;
	string input_file;
	flt concentrations;
	
	flt dpd_density;
//	vector<int> place;
	vect3 box_size;
	
	int steps_1;
	int steps_2;
	int steps_3;
	flt timestep;
	
	vector<vector<flt> > pair_coeff;
	
	int output_freq;
	
	vector<vector<flt> > bond_template;
	int steps_ch;
	
//	vector<flt> chem_tri;
//	vector<int> chem_para;
	
	vector<string> bead_type;
	int restart;
	
	vector<vector<vector<flt> > > angle_templateAngle;
	vector<vector<vector<flt> > > angle_templateForce;
};

Conf readDpdConf(string fname);
vector<string> parseLine(string line);

#endif //LPD_CONF_H
