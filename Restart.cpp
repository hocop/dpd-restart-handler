//
// Created by rusl on 3/21/17.
//

#include "Restart.h"
#include <iomanip>

using namespace std;

void distPBC_c(flt &dx, flt &dy, flt &dz, flt const &gridx, flt const &gridy, flt const &gridz) {
	dx = (fabs(dx) < gridx / 2.0) ? dx : (dx > 0) ? dx - gridx : dx + gridx;
	dy = (fabs(dy) < gridy / 2.0) ? dy : (dy > 0) ? dy - gridy : dy + gridy;
	dz = (fabs(dz) < gridz / 2.0) ? dz : (dz > 0) ? dz - gridz : dz + gridz;
}
void absPBC_c_while(flt &x, flt &y, flt &z, flt const &gridx, flt const &gridy, flt const &gridz) {
	while (x < 0.f)
		x += gridx;
	while (x > gridx)
		x -= gridx;
	while (y < 0.f)
		y += gridy;
	while (y > gridy)
		y -= gridy;
	while (z < 0.f)
		z += gridz;
	while (z > gridz)
		z -= gridz;
}
void absPBC_c_if(flt &x, flt &y, flt &z, flt const &gridx, flt const &gridy, flt const &gridz) {
	if(x < 0.f)
		x += gridx;
	else if(x > gridx)
		x -= gridx;
	if(y < 0.f)
		y += gridy;
	else if(y > gridy)
		y -= gridy;
	if(z< 0.f)
		z += gridz;
	else if(z > gridz)
		z -= gridz;
}

void Restart::distPBC(vect3 &x) const {
	distPBC_c(x[0], x[1], x[2], box[0], box[1], box[2]);
}

void Restart::absPBC_while(vect3 &x) const {
	absPBC_c_while(x[0], x[1], x[2], box[0], box[1], box[2]);
}
void Restart::absPBC_if(vect3 &x) const {
	absPBC_c_if(x[0], x[1], x[2], box[0], box[1], box[2]);
}

Restart::Restart() {
	natoms = nbonds = nchains = nangles = nnodes = 0;
	x = vector<vect3>(0);
	valency = vector<int>(0);
	type = vector<int>(0);
	flag = vector<int>(0);
	a1 = vector<int>(0);
	a2 = vector<int>(0);
	
	chain = vector<vector<int> >(0);
	center = vector<int>(0);
	shell = vector<vector<int> >(0);
	bond = vector<vector<int> >(0);
}

Restart::Restart(Conf const &conf) : Restart() {
	box = conf.box_size;
	density = conf.dpd_density;
}

Restart::Restart(Restart const &r) {
	natoms = r.natoms;
	nbonds = r.nbonds;
	nchains = r.nchains;
	nangles = r.nangles;
	nnodes = r.nnodes;
	
	density = r.density;
	
	x = vector<vect3>(natoms);
	valency = vector<int>(natoms);
	type = vector<int>(natoms);
	flag = vector<int>(natoms);
	for (int i = 0; i < natoms; ++i) {
		x[i] = r.x[i];
		valency[i] = r.valency[i];
		type[i] = r.type[i];
		flag[i] = r.flag[i];
	}
	a1 = vector<int>(nbonds);
	a2 = vector<int>(nbonds);
	for (int i = 0; i < nbonds; ++i) {
		a1[i] = r.a1[i];
		a2[i] = r.a2[i];
	}
	box = r.box;
	
	chain = vector<vector<int> >(nchains);
	center = vector<int>(nnodes);
	shell = vector<vector<int> >(nnodes);
	bond = vector<vector<int> >(nnodes);
	for (int i = 0; i < nchains; ++i) {
		chain[i] = vector<int>(r.chain[i].size());
		for (int j = 0; j < r.chain[i].size(); ++j)
			chain[i][j] = r.chain[i][j];
	}
	for(int i = 0; i < nnodes; ++i) {
		center[i] = r.center[i];
		shell[i] = vector<int>(r.shell[i].size());
		for (int j = 0; j < r.shell[i].size(); ++j)
			shell[i][j] = r.shell[i][j];
		bond[i] = vector<int>(r.bond[i].size());
		for (int j = 0; j < r.bond[i].size(); ++j)
			bond[i][j] = r.bond[i][j];
	}
}

void Restart::add(Restart const &r, vect3 r0, flt a, flt b, flt c, flt scale) {
	for(int i = 0; i < r.natoms; ++i) {
		vect3 xn = r0 + scale * rotate(r.x[i], a, b, c);
		absPBC_while(xn);
		x.push_back(xn);
		valency.push_back(r.valency[i]);
		type.push_back(r.type[i]);
		flag.push_back(r.flag[i]);
	}
	for(int i = 0; i < r.nbonds; ++i) {
		a1.push_back(r.a1[i] + natoms);
		a2.push_back(r.a2[i] + natoms);
	}
	for(int i = 0; i < r.nchains; ++i) {
		vector<int> ch(r.chain[i].size());
		for(int j = 0; j < ch.size(); ++j)
			ch[j] = r.chain[i][j] + natoms;
		chain.push_back(ch);
	}
	for(int i = 0; i < r.nnodes; ++i) {
		center.push_back(r.center[i] + natoms);
		vector<int> sh(r.shell[i].size());
		for(int j = 0; j < sh.size(); ++j)
			sh[j] += natoms;
		shell.push_back(sh);
		vector<int> b(r.bond[i].size());
		for(int j = 0; j < b.size(); ++j)
			b[j] += nbonds;
		bond.push_back(b);
	}
	
	natoms += r.natoms;
	nbonds += r.nbonds;
	nangles += r.nangles;
	nchains += r.nchains;
	nnodes += r.nnodes;
}

void topology(int natoms, int nbonds, vector<int> const &a1, vector<int> const &a2, vector<vector<int> > &chains_Id,
			  vector<int> &center, vector<vector<int> > &shell, vector<vector<int> > &bond, int &nangles)
{
	//!topology analysis
	std::vector<int> bpa(natoms); // bonds per atom
	for (int i = 0; i < natoms; ++i)
		bpa[i] = 0;
	for (int i = 0; i < nbonds; ++i) {
		++bpa[a1[i]];
		++bpa[a2[i]];
	}
	std::vector<int> whatAtomIs(natoms);
	/**
	 * -1: not stated
	 * 0: free particle
	 * 1: node
	 * 2: chain part
	 * 20: also chain part (used to mark atoms that are already added to chain)
	 */
	for (int i = 0; i < natoms; ++i)
		whatAtomIs[i] = -1;
	for (int i = 0; i < natoms; ++i) {
		whatAtomIs[i] = (bpa[i] == 2) ? 2 : (bpa[i] == 0) ? 0 : (bpa[i] != 1) ? 1 : 0;
		if (bpa[i] == 1) { // pairs are viewed separately
			int otherId = -1;
			for (int j = 0; otherId == -1; ++j) // (j < nbonds) automatically
				otherId = (a1[j] == i) ? a2[j] : (a2[j] == i) ? a1[j] : otherId;
			if (bpa[otherId] == 1 && i > otherId)
				whatAtomIs[i] = 1; // make one atom of two a node
		}
	}
	
	std::vector<bool> bondCounted(nbonds); // check if this bond is already being computed by chain or node
	for (int i = 0; i < nbonds; ++i)
		bondCounted[i] = false;
	
	//! construct chains
	for (int i = 0; i < natoms; ++i) {
		if (whatAtomIs[i] == 2) {
			std::vector<int> chain; // chain to be added to "chains"
			
			// collect info on this chain link
			int left = -1, right = -1;
			for (int j = 0; j < nbonds; ++j) {
				int li = (a1[j] == i) ? a2[j] : (a2[j] == i) ? a1[j] : -1;
				int &wl = (left == -1) ? left : right;
				wl = (li != -1) ? li : wl;
			}
			
			// make two parts to connect them later
			std::vector<int> chain1, chain2;
			
			whatAtomIs[i] = 20;
			for (int dir = 0; dir < 2; ++dir) {
				// go one direction of two
				int ci = (dir == 0) ? left : right;
				int lasti = i;
				while (whatAtomIs[ci] == 2) {
					if (dir == 0)
						chain1.push_back(ci);
					else
						chain2.push_back(ci);
					//go to next atom
					int ci_ = ci;
					for (int j = 0; j < natoms; ++j) {
						if (a1[j] == ci && a2[j] != lasti) {
							ci = a2[j];
							bondCounted[j] = true; // included in chain
							break;
						}
						if (a2[j] == ci && a1[j] != lasti) {
							ci = a1[j];
							bondCounted[j] = true; // included in chain
							break;
						}
					}
					lasti = ci_;
					whatAtomIs[lasti] = 20;
				}
				// add ends
				if (dir == 0)
					chain1.push_back(ci);
				else
					chain2.push_back(ci);
			}
			
			//connect the two parts of chain
			std::reverse(chain1.begin(), chain1.end());
			chain1.push_back(i);
			chain = chain1;
			if (chain[0] != i) // chain is not cyclic
				chain.insert(chain.end(), chain2.begin(), chain2.end());
			else // chain is cyclic
				whatAtomIs[i] = 1;
			
			//finally
			chains_Id.push_back(chain);
		}
	}
	
	
	// nodes
	int nnodes = 0;
	for (int i = 0; i < natoms; ++i)
		if (whatAtomIs[i] == 1) // is node
			++nnodes;
	center.resize(nnodes);
	shell.resize(nnodes);
	bond.resize(nnodes);
	
	// construct nodes
	int ia = 0;
	for (int i = 0; i < nnodes; ++i) {
		// find next node center atom
		while (whatAtomIs[ia] != 1)
			++ia;
		center[i] = ia;
		
		// make shell
		int shSize = 0;
		for (int j = 0; j < natoms; ++j)
			if (a1[j] == ia || a2[j] == ia)
				++shSize;
		shell[i].resize(shSize);
		int iS = 0;
		for (int j = 0; j < natoms; ++j) {
			shell[i][iS] = (a1[j] == ia) ? a2[j] : (a2[j] == ia) ? a1[j] : -1;
			iS += (shell[i][iS] == -1) ? 0 : 1;
		}
		++ia;
	}
	// tell nodes what bonds to compute
	for (int i = 0; i < nnodes; ++i) {
		int centId = center[i];
		for (int j = 0; j < shell[i].size(); ++j) {
			int atomId = shell[i][j];
			int what = whatAtomIs[atomId];
			
			if (what == 0) // free atom
				bond[i].push_back(atomId);
			
			if (what == 1) { // another node
				int bondId = -1; // find bond that connects these nodes
				for (int k = 0; k < nbonds; ++k)
					if ((a1[k] == atomId && a2[k] == centId) || (a1[k] == centId && a2[k] == atomId)) {
						bondId = k;
						break;
					}
				if (bondCounted[bondId])
					continue;
				bondCounted[bondId] = true;
				// find node that is that shell atom
				int nodeId = -1;
				for (int k = 0; k < nnodes; ++k)
					if (center[k] == atomId) {
						nodeId = k;
						break;
					}
				// my tricky condition:
				if ((nodeId + nnodes - i) % nnodes > nnodes / 2)
					continue; // for uniform distribution of computations over nodes
				if ((shell[i].size() <= shell[nodeId].size()) || (bond[i].size() < bond[nodeId].size()))
					bond[i].push_back(atomId); // this atom will comp bond
				else
					bond[nodeId].push_back(centId); // another atom will compute this bond
			}
			
			// if what == 20 then nothing - chain computes itself
		}
	}
	nangles = 0;
	for (int i = 0; i < chains_Id.size(); ++i)
		nangles += (chains_Id[i].size() - 2 > 0) ? chains_Id[i].size() - 2 : 0;
	for (int i = 0; i < center.size(); ++i)
		nangles += shell[i].size() * (shell[i].size() - 1);
}

void topology(Restart &r) {
	topology(r.natoms, r.nbonds, r.a1, r.a2, r.chain, r.center, r.shell, r.bond, r.nangles);
	r.nchains = (int) r.chain.size();
	r.nnodes = (int) r.center.size();
}

void topologyReverse(Restart &r) { // unchecked
	r.a1.clear();
	r.a2.clear();
	for(int i = 0; i < r.nchains; ++i)
		for(int j = 0; j < r.chain[i].size() - 1; ++j) {
			r.a1.push_back(r.chain[i][j]);
			r.a2.push_back(r.chain[i][j+1]);
		}
	for(int i = 0; i < r.nnodes; ++i)
		for(int j = 0; j < r.shell[i].size(); ++j) {
			r.a1.push_back(r.center[i]);
			r.a2.push_back(r.shell[i][j]);
		}
	r.nbonds = (int) r.a1.size();
}
Restart importGav(string fname, int withflag, bool withTopology) {
	ifstream f;
	f.open(fname.c_str());
	int natoms;
	double gX, gY, gZ;
	flt ppc; //particle per cell
	f >> natoms >> ppc;
	f >> gX >> gY >> gZ;
	
	// list of things that are decremented:
	// index of atom
	// type of atom
	
	vector<int> valency(natoms);
	vector<int> type(natoms);
	vector<int> flag(natoms);
	vector<vect3> x(natoms);
	for (int i = 0; i < natoms; ++i) {
		int n, val, fl;
		f >> n;
		--n;
		int type_tmp;
		f >> val >> type_tmp;
		if(withflag == 1)
			f >> fl;
		f >> x[n][0] >> x[n][1] >> x[n][2];
		type[n] = type_tmp - 1;
		valency[n] = val;
		flag[n] = fl;
	}
	
	string str_bonds;
	int nbonds;
	f >> str_bonds >> nbonds;
	vector<int> a1(nbonds), a2(nbonds);
	for (int i = 0; i < nbonds; ++i) {
		int aa1, aa2;
		f >> aa1 >> aa2;
		a1[i] = aa1 - 1;
		a2[i] = aa2 - 1;
	}
	f.close();
	
	std::vector<std::vector<int> > chains_Id;
	vector<int> center;
	vector<vector<int> > shell;
	std::vector<std::vector<int> > bond;
	int nangles;
	if (withTopology)
		topology(natoms, nbonds, a1, a2, chains_Id, center, shell, bond, nangles);
	
	Restart rst;
	rst.natoms = (int) x.size();
	rst.nbonds = (int) a1.size();
	rst.nangles = nangles;
	rst.nchains = chains_Id.size();
	rst.nnodes = center.size();
	rst.density = ppc;
	
	rst.x = x;
	rst.valency = valency;
	rst.type = type;
	rst.flag = flag;
	rst.a1 = a1;
	rst.a2 = a2;
	rst.box = vect3(gX, gY, gZ);
	
	rst.chain = chains_Id;
	rst.center = center;
	rst.shell = shell;
	rst.bond = bond;
	return rst;
}

Restart importMol(string fname, Conf conf) {
	ifstream f;
	f.open(fname.c_str());
	Restart r(conf);
	
	string none;
	f >> none;
	f >> none;
	f >> r.natoms >> r.nbonds;
	for(int j = 0; j < 9; ++j)
		f >> none;
	r.x = vector<vect3>(r.natoms);
	r.valency = vector<int>(r.natoms);
	r.type = vector<int>(r.natoms);
	r.flag = vector<int>(r.natoms);
	r.a1 = vector<int>(r.nbonds);
	r.a2 = vector<int>(r.nbonds);
	for(int i = 0; i < r.natoms; ++i) {
		string t;
		f >> r.x[i][0] >> r.x[i][1] >> r.x[i][2] >> t;
		int ti = -1;
		for(int j = 0; j < conf.bead_type.size(); ++j)
			if(conf.bead_type[j] == t){
				ti = j;
				break;
			}
		r.type[i] = ti;
		for(int j = 0; j < 12; ++j)
			f >> none;
	}
	for(int i = 0; i < r.nbonds; ++i){
		f >> r.a1[i] >> r.a2[i];
		--r.a1[i];
		--r.a2[i];
		for(int j = 0; j < 5; ++j)
			f >> none;
	}
	f.close();
	topology(r);
	return r;
}

#define BASE 223

unsigned char intToUC(int a) { // a from 0 to (BASE-1)
	return (unsigned char) (a + (256-BASE));
}
int ucToInt(unsigned char a) {
	return (int) a - (256-BASE);
}

string intToSmall(int a, int digits) {
	vector<int> d(digits);
	string out = "";
	for(int i = 0; i < digits; ++i, a /= BASE)
		out = string(1, intToUC(a % BASE)) + out;
	return out;
}

int smallToInt(string s, int digits) {
	int out = 0, ex = 1;
	for(int i = 0; i < digits; ++i, ex *= BASE)
		out += ex * ucToInt(s[digits-1-i]);
	return out;
}

string fltToSmall(flt a, int digits) {
	return intToSmall((int)a, digits-1) + intToSmall((int)round((a-(int)a)*(BASE-1)), 1);
}

flt smallToFlt(string s, int digits) {
	string a("");
	for(int i = 0; i < s.size()-1; ++i)
		a+=s[i];
	string b = string("") + s[s.size()-1];
	return smallToInt(a, digits-1) + smallToInt(b, 1) / (flt)(BASE-1);
}

vector<int> lineToSetInt(string s, int digits) {
	int n = (int) s.size() / digits;
	vector<int> res;
	for(int i = 0; i < n; ++i) {
		string ss("");
		for(int j = 0; j < digits; ++j)
			ss += string(1, s[i * digits + j]);
		res.push_back(smallToInt(ss, digits));
	}
	return res;
}

vector<flt> lineToSetFlt(string s, int digits) {
	int n = (int) s.size() / digits;
	vector<flt> res;
	for(int i = 0; i < n; ++i) {
		string ss("");
		for(int j = 0; j < digits; ++j)
			ss += string(1, s[i * digits + j]);
		res.push_back(smallToFlt(ss, digits));
	}
	return res;
}

Restart importSmall(string fname) {
	ifstream f;
	f.open(fname.c_str());
	Restart r;
	
	string none;
	f >> none >> r.density;
	f >> none >> r.box.x >> r.box.y >> r.box.z;
	f >> none >> r.natoms;
	f >> none >> r.nbonds;
	f >> none >> r.nchains;
	f >> none >> r.nnodes;
	
	int digits = (r.box.x > BASE-1 || r.box.y > BASE-1 || r.box.z > BASE-1)? 3 : 2;
	int intDigits = (int) (log(r.natoms) / log(BASE)+.001f) + 1;
	// read all atoms info
	r.valency = vector<int>(r.natoms);
	r.type = vector<int>(r.natoms);
	r.flag = vector<int>(r.natoms);
	string line;
	f >> line;
	vector<int> all = lineToSetInt(line, 1);
	for(int i = 0; i < r.natoms; ++i) {
		r.valency[i] = all[3*i];
		r.type[i] = all[3*i+1];
		r.flag[i] = all[3*i+2];
	}
	// read atoms cords
	f >> line;
	vector<flt> set = lineToSetFlt(line, digits);
	r.x = vector<vect3>(r.natoms);
	for(int i = 0; i < r.natoms; ++i)
		r.x[i] = vect3(set[i*3], set[i*3+1], set[i*3+2]);
	// read chains
	f >> line;
	all = lineToSetInt(line, intDigits);
	for(int i = 0; i < all.size();) {
		int size = all[i];
		vector<int> chain(size);
		for(int j = 0; j < size; ++j)
			chain[j] = all[i+j+1];
		r.chain.push_back(chain);
		i += size + 1;
	}
	// read nodes
	if(r.nnodes != 0) {
		f >> line;
		all = lineToSetInt(line, intDigits);
		for (int i = 0; i < all.size();) {
			int size = all[i];
			r.center.push_back(all[i + 1]);
			vector<int> shell(size);
			for (int j = 0; j < size; ++j)
				shell[j] = all[i + j + 2];
			r.shell.push_back(shell);
			i += size + 2;
		}
	}
	f.close();
	topologyReverse(r);
	return r;
}

Restart importNormal(string fname) {
	ifstream f;
	f.open(fname.c_str());
	Restart r;
	
	string none;
	f >> none >> r.density;
	f >> none >> r.box.x >> r.box.y >> r.box.z;
	f >> none >> r.natoms;
	f >> none >> r.nbonds;
	f >> none >> r.nchains;
	f >> none >> r.nnodes;
	
	// read all atoms info and cords
	f >> none; // word "Atoms"
	r.valency = vector<int>(r.natoms);
	r.type = vector<int>(r.natoms);
	r.flag = vector<int>(r.natoms);
	r.x = vector<vect3>(r.natoms);
	for(int i = 0; i < r.natoms; ++i)
		f >> r.valency[i] >> r.type[i] >> r.flag[i] >> r.x[i].x >> r.x[i].y >> r.x[i].z;
	// read chains
	f >> none; // word "Chains"
	for(int i = 0; i < r.nchains; ++i) {
		int size;
		f >> size;
		vector<int> chain(size);
		for (int j = 0; j < size; ++j)
			f >> chain[j];
		r.chain.push_back(chain);
	}
	// read nodes
	f >> none; // word "Nodes"
	for(int i = 0; i < r.nnodes; ++i) {
		int center, size;
		f >> center >> size;
		r.center.push_back(center);
		vector<int> shell(size);
		for(int j = 0; j < size; ++j)
			f >> shell[j];
		r.shell.push_back(shell);
	}
	f.close();
	topologyReverse(r);
	return r;
}

void writeval(int val, std::ostream &out, int tabsize) { // this is not mine
	int vc(val), dval(1);
	while (vc >= 10 || vc <= -10) {
		dval++;
		vc /= 10;
	}
	//if(val<0) dval++;
	for (int q = dval; q < tabsize; q++) out << " ";
	out << val;
}

void writeval(double val, int precision, std::ostream &out, int tabsize) { // this is not mine
	double vc(val);
	int dval(1);
	while (vc >= 10.00 || vc <= -10.00) {
		dval++;
		vc /= 10.00;
	}
	//if(val<0.00) dval++;
	for (int q = dval; q < tabsize; q++) out << " ";
//	out<<v3_general::dtos_sp(val, precision);
	out << setprecision(precision);
	out << val;
}

void exportGav(Restart &rst, string fname, int withflag) {
	ofstream out(fname.c_str());
	writeval(rst.natoms, out, TABSIZE);
	writeval(rst.density, 4, out, 2);
	out << "\n";
	writeval(rst.box.x, 12, out, 6);
	writeval(rst.box.y, 12, out, 6);
	writeval(rst.box.z, 12, out, 6);
	for (int q = 0; q < rst.natoms; q++) {
		out << "\n";
		writeval(q + 1, out, TABSIZE);
		writeval(rst.valency[q], out, 4);
		writeval(rst.type[q] + 1, out, 4);
		if(withflag == 1)
			writeval(rst.flag[q], out, 4);
		writeval(rst.x[q][0], 6, out, 6);
		writeval(rst.x[q][1], 6, out, 6);
		writeval(rst.x[q][2], 6, out, 6);
	}
	out << "\n bonds:";
	writeval(rst.nbonds, out, 11);
	for (int q = 0; q < rst.nbonds; q++) {
		out << "\n";
		writeval(rst.a1[q] + 1, out, TABSIZE);
		writeval(rst.a2[q] + 1, out, TABSIZE);
	}
	out << "\n angles:";
	writeval(rst.nangles, out, 10);
	for(int q = 0; q < rst.nchains; ++q)
		for(int i = 0; i < rst.chain[q].size()-2; ++i) {
			out << "\n";
			writeval(rst.chain[q][i] + 1, out, TABSIZE);
			writeval(rst.chain[q][i + 1] + 1, out, TABSIZE);
			writeval(rst.chain[q][i + 2] + 1, out, TABSIZE);
		}
	for(int q = 0; q < rst.nnodes; ++q)
		for(int i = 0; i < rst.shell[q].size(); ++i)
			for(int j = i + 1; j < rst.shell[q].size(); ++j){
				out << "\n";
				writeval(rst.shell[q][i] + 1, out, TABSIZE);
				writeval(rst.center[q] + 1, out, TABSIZE);
				writeval(rst.shell[q][j] + 1, out, TABSIZE);
			}
	out.close();
}

void exportSmall(Restart &r, string fname) {
	ofstream f(fname.c_str());
	f << "dpd_density\t" << r.density << endl;
	f << "box_size\t" << r.box.x << '\t' << r.box.y << '\t' << r.box.z << endl;
	f << "atoms_number\t" << r.natoms << endl;
	f << "bonds_number\t" << r.nbonds << endl;
	f << "chains_number\t" << r.nchains << endl;
	f << "nodes_number\t" << r.nnodes << endl;
	int digits = (r.box.x > BASE-1 || r.box.y > BASE-1 || r.box.z > BASE-1)? 3 : 2;
	int intDigits = (int) (log(r.natoms) / log(BASE)+.001f) + 1;
	// write all atoms info
	// todo: assume these numbers small
	for(int i = 0; i < r.natoms; ++i)
		f << intToSmall(r.valency[i], 1) << intToSmall(r.type[i], 1) << intToSmall(r.flag[i], 1);
	f << endl;
	// write all atoms cords
	for(int i = 0; i < r.natoms; ++i)
		f << fltToSmall(r.x[i][0], digits) << fltToSmall(r.x[i][1], digits) << fltToSmall(r.x[i][2], digits);
	f << endl;
	// write chains
	for(int i = 0; i < r.nchains; ++i) {
		f << intToSmall(r.chain[i].size(), intDigits);
		for (int j = 0; j < r.chain[i].size(); ++j)
			f << intToSmall(r.chain[i][j], intDigits);
	}
	f << endl;
	// write all nodes
	for(int i = 0; i < r.nnodes; ++i) {
		f << intToSmall(r.center[i], intDigits);
		f << intToSmall(r.shell[i].size(), intDigits);
		for(int j = 0; j < r.shell[i].size(); ++j)
			f << intToSmall(r.shell[i][j], intDigits);
	}
	f.close();
}

void exportNormal(Restart &r, string fname) {
	ofstream f(fname.c_str());
	f << "dpd_density\t" << r.density << endl;
	f << "box_size\t" << r.box.x << '\t' << r.box.y << '\t' << r.box.z << endl;
	f << "atoms_number\t" << r.natoms << endl;
	f << "bonds_number\t" << r.nbonds << endl;
	f << "chains_number\t" << r.nchains << endl;
	f << "nodes_number\t" << r.nnodes << endl;
	// write all atoms info and cords
	f << "Atoms:\n";
	for(int i = 0; i < r.natoms; ++i)
		f << r.valency[i] << ' ' << r.type[i] << ' ' << r.flag[i] << '\t' << r.x[i] << '\n';
	f << endl;
	// write chains
	f << "Chains:\n";
	for(int i = 0; i < r.nchains; ++i) {
		f << r.chain[i].size() << '\n';
		for (int j = 0; j < r.chain[i].size(); ++j)
			f << ' ' << r.chain[i][j];
		f << '\n';
	}
	// write all nodes
	f << "Nodes:\n";
	for(int i = 0; i < r.nnodes; ++i) {
		f << r.center[i] << '\t';
		f << r.shell[i].size() << '\n';
		for(int j = 0; j < r.shell[i].size(); ++j)
			f << r.shell[i][j] << ' ';
		f << '\n';
	}
	f.close();
}


















