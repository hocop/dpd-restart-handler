//
// Created by rusl on 3/21/17.
//

#include <fstream>
#include "Conf.h"

vector<string> parseLine(string line) {
	vector<string> parsed;
	bool split = true;
	string temp = "";
	
	for( int i = 0; i < line.length(); ++i ) {
		char c = line[i];
		if(c == '#') // comment
			break;
		if( c == ' ' || c == '\t' ) {
			if( !split ) {
				parsed.push_back(temp);
				temp = "";
			}
			split = true;
		} else {
			split = false;
			temp += c;
		}
	}
	if( !split )
		parsed.push_back(temp);
	
	return parsed;
}

Conf readDpdConf(string fname) {
	ifstream f(fname.c_str());
	if( !f.is_open() ) {
		throw "Error: could not open conf file\n";
	}
	
	string line;
	vector<string> parsed;
	
	Conf conf;
	
	int TYPES_N = -1;
	bool typesread = false;
	
	getline(f, line);
	parsed = parseLine( line );
	while(1) {
		getline(f, line);
		parsed = parseLine( line );
		if(parsed.size() == 0)
			continue;
		if(parsed[0] == "end")
			break;
		
		if(parsed[0] == "run_type")
			conf.run_type = (parsed[1] == "re")? re : co;
		if(parsed[0] == "species_number") {
			conf.species_number = atoi(parsed[1].c_str());
			TYPES_N = 4; //todo: HOW TO READ?
			conf.bead_type = vector<string>(TYPES_N);
			conf.pair_coeff = vector<vector<flt> >(TYPES_N);
			conf.bond_template = vector<vector<flt> >(TYPES_N);
			conf.angle_templateAngle = vector<vector<vector<flt> > >(TYPES_N);
			conf.angle_templateForce = vector<vector<vector<flt> > >(TYPES_N);
			for(int j = 0; j < TYPES_N; ++j) {
				conf.pair_coeff[j] = vector<flt>(TYPES_N);
				conf.bond_template[j] = vector<flt>(TYPES_N);
				conf.angle_templateAngle[j] = vector<vector<flt> >(TYPES_N);
				conf.angle_templateForce[j] = vector<vector<flt> >(TYPES_N);
				for(int k = 0; k < TYPES_N; ++k) {
					conf.angle_templateAngle[j][k] = vector<flt>(TYPES_N);
					conf.angle_templateForce[j][k] = vector<flt>(TYPES_N);
				}
			}
		}
		if(parsed[0] == "input_file")
			conf.input_file = parsed[1];
		if(parsed[0] == "concentrations")
			conf.concentrations = (flt) atof(parsed[1].c_str());
		
		if(parsed[0] == "dpd_density")
			conf.dpd_density = (flt) atof(parsed[1].c_str());
		if(parsed[0] == "box_size")
			conf.box_size = vect3(atoi(parsed[1].c_str()), atoi(parsed[2].c_str()), atoi(parsed[3].c_str()));
		
		if(parsed[0] == "steps_1")
			conf.steps_1 = atoi(parsed[1].c_str());
		if(parsed[0] == "steps_2")
			conf.steps_2 = atoi(parsed[1].c_str());
		if(parsed[0] == "steps_3")
			conf.steps_3 = atoi(parsed[1].c_str());
		if(parsed[0] == "timestep" )
			conf.timestep = (flt) atof(parsed[1].c_str());
		
		if(parsed[0] == "pair_coeff"){
			if(TYPES_N < 0)
				throw 1;
			// todo: fill
		}
		
		if(parsed[0] == "output_freq")
			conf.output_freq = atoi(parsed[1].c_str());
		
		if(parsed[0] == "bond_template"){
			if(TYPES_N < 0)
				throw 1;
			// todo: fill
		}
		if(parsed[0] == "steps_ch")
			conf.steps_ch = atoi(parsed[1].c_str());
		
		if(parsed[0] == "bead_type")
			if(TYPES_N < 0)
				throw 1;
			else
				conf.bead_type[atoi(parsed[1].c_str())-1] = parsed[2];
		if(parsed[0] == "r")
			conf.restart = atoi(parsed[1].c_str());
		if(parsed[0] == "angle_template"){
			if(TYPES_N < 0)
				throw 1;
			// todo: fill
		}
	}
	
	/* information has been read, now init dpdconf */
	return conf;
}
