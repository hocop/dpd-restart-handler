//
// Created by rusl on 4/11/17.
//

#include "Analyser.h"

bool contains(vector<int> &v, int e) {
	for(int i = 0; i < v.size(); ++i)
		if(v[i] == e)
			return true;
	return false;
}

Restart selectByType(Restart const &r, vector<int> &type) {
	Restart k(r);
	k.x.clear();
	k.valency.clear();
	k.type.clear();
	k.flag.clear();
	
	vector<int> shift(r.natoms);
	for(int i = 0; i < r.natoms; ++i)
		shift[i] = 0;
	for(int i = 0; i < r.natoms; ++i)
		if(contains(type, r.type[i])) {
			k.x.push_back(r.x[i]);
			k.valency.push_back(r.valency[i]);
			k.type.push_back(r.type[i]);
			k.flag.push_back(r.flag[i]);
		} else
			for(int j = i+1; j < r.natoms; ++j)
				--shift[j];
	
	k.a1.clear();
	k.a2.clear();
	for(int i = 0; i < r.nbonds; ++i)
		if (contains(type, r.type[r.a1[i]]) && contains(type, r.type[r.a2[i]])) {
			k.a1.push_back(r.a1[i] + shift[r.a1[i]]);
			k.a2.push_back(r.a2[i] + shift[r.a2[i]]);
		}
	
	k.natoms = (int) k.x.size();
	k.nbonds = (int) k.a1.size();
	
	k.chain.clear();
	k.shell.clear();
	k.center.clear();
	k.bond.clear();
	topology(k);
	k.nchains = k.chain.size();
	k.nnodes = k.center.size();
	
	int nangles = 0;
	for (int i = 0; i < k.chain.size(); ++i)
		nangles += (k.chain[i].size() - 2 > 0) ? k.chain[i].size() - 2 : 0;
	for (int i = 0; i < k.center.size(); ++i)
		nangles += k.shell[i].size() * (k.shell[i].size() - 1);
	
	k.nangles = nangles;
	
	return k;
}

vector<flt> countNematicOrderBonds(Restart const &r, int rmin, int rmax) {
	vector<vect3> dir(r.nbonds); // direction of every bond
	for(int i = 0; i < r.nbonds; ++i) {
		dir[i] = r.x[r.a2[i]] - r.x[r.a1[i]];
		r.distPBC(dir[i]);
		dir[i] /= sqrt(dir[i] * dir[i]);
	}
	vector<flt> Sr(rmax - rmin);
	vector<flt> nCounted(rmax - rmin);
	for (int ri = rmin; ri < rmax; ++ri) {
		Sr[ri - rmin] = 0.0;
		nCounted[ri - rmin] = 0;
	}
	int accuracy = r.nbonds / 500;
	for (int i = 0; i < r.nbonds; i += accuracy) {
		vector<vect3> mdir(rmax - rmin); // medium direction {x, y, z}
		for (int j = 0; j < r.nbonds; ++j) {
			vect3 dx = r.x[r.a1[i]] - r.x[r.a1[j]];
			r.distPBC(dx);
			flt rsq = dx * dx;
			for (int ri = rmax - 1; ri >= rmin; --ri) {
				flt rsqOk = ri * ri;
				if (rsq < rsqOk) {
					flt mult = (mdir[ri - rmin] * dir[j] > 0) ? 1.f : -1.f;
					mdir[ri - rmin] += mult * dir[j];
				} else
					break;
			}
		}
		for (int ri = rmin; ri < rmax; ++ri) {
			int i = ri - rmin;
			mdir[i] /= sqrt(mdir[i] * mdir[i]);
		}
		// directors computed, compute S(r)
		for (int j = 0; j <r.nbonds; ++j) {
			vect3 dx = r.x[r.a1[i]] - r.x[r.a1[j]];
			r.distPBC(dx);
			flt rsq = dx * dx;
			for (int ri = rmax - 1; ri >= rmin; --ri) {
				flt rsqOk = ri * ri;
				if (rsq < rsqOk) {
					flt cocos = dir[i] * mdir[ri - rmin];
					cocos *= cocos;
					Sr[ri - rmin] += cocos;
					++nCounted[ri - rmin];
				} else
					break;
			}
		}
	}
	
	for (int ri = rmin; ri < rmax; ++ri) {
		Sr[ri - rmin] /= nCounted[ri - rmin];
		Sr[ri - rmin] = .5 * (3 * Sr[ri - rmin] - 1);
	}
	
	return Sr;
}

vector<flt> countNematicOrderChains(Restart const &r, int rmin, int rmax) {
	++rmax;
	vector<vect3> dir(r.nchains); // direction of every chain
	for(int i = 0; i < r.nchains; ++i) {
		dir[i] = r.x[r.chain[i][r.chain[i].size()-1]] - r.x[r.chain[i][0]];
		r.distPBC(dir[i]);
		dir[i] /= dir[i].mod();
	}
	
	vector<flt> Sr(rmax - rmin);
	vector<flt> nCounted(rmax - rmin);
	for(int ri = rmin; ri < rmax; ++ri) {
		Sr[ri - rmin] = 0.f;
		nCounted[ri - rmin] = 0;
	}
	int accuracy = 1 + r.nchains / 500;
	for(int i = 0; i < r.nchains; i += accuracy) {
		vector<vect3> mdir(rmax - rmin);
		for (int j = 0; j < r.nchains; j += 1) {
			vect3 dx = centerOfMassOfChain(r, i) - centerOfMassOfChain(r, j);
			r.distPBC(dx);
			flt rsq = dx * dx;
			for (int ri = rmax - 1; ri >= rmin; --ri) {
				flt rsqOk = ri * ri;
				if (rsq < rsqOk) {
					mdir[ri-rmin] += (dir[j] * mdir[ri-rmin] > 0)? dir[j] : -dir[j];
				} else
					break;
			}
		}
		for(int ri = rmin; ri < rmax; ++ri)
			mdir[ri - rmin] /= mdir[ri - rmin].mod();
		for (int j = 0; j < r.nchains; j += 1) {
			if(i == j)
				continue;
			vect3 dx = centerOfMassOfChain(r, i) - centerOfMassOfChain(r, j);
			r.distPBC(dx);
			flt rsq = dx * dx;
			for (int ri = rmax - 1; ri >= rmin; --ri) {
				flt cocos = mdir[ri-rmin] * dir[j];
				cocos *= cocos;
				flt rsqOk = ri * ri;
				if (rsq < rsqOk) {
					Sr[ri - rmin] += cocos;
					++nCounted[ri - rmin];
				} else
					break;
			}
		}
	}
	
	for(int ri = rmin; ri < rmax; ++ri) {
		Sr[ri - rmin] /= nCounted[ri - rmin];
		Sr[ri - rmin] = .5f * (3.f*Sr[ri - rmin] - 1.f);
	}
	
	return Sr;
}

flt countNematicOrderChainsMaxR(Restart const &r) {
	int sz = (int) r.box.normMax(); // size of box
	return countNematicOrderChains(r, sz, sz)[0];
}

void meanAndStdev(vector<flt> value, flt &mean, flt &stdev) {
	int n = (int) value.size();
	mean = stdev = 0.f;
	for(int i = 0; i < n; ++i)
		mean += value[i];
	mean /= n;
	for(int i = 0; i < n; ++i) {
		flt d = value[i] - mean;
		stdev += d*d;
	}
	stdev /= n - 1;
	stdev = sqrt(stdev);
}

void countEEDist(Restart const &r, flt &mean, flt &stdev) {
	vector<flt> dist;
	for(int i = 0; i < r.chain.size(); ++i) {
		vect3 dx = r.x[r.chain[i][0]] -
					r.x[r.chain[i][r.chain[i].size() - 1]];
		r.distPBC(dx);
		dist.push_back(dx.mod());
	}
	meanAndStdev(dist, mean, stdev);
}

vect3 countDirectorChains(Restart const &r) {
	vect3 dir(0, 0, 0);
	for(int i = 0; i < r.chain.size(); ++i) {
		vect3 d = (r.x[r.chain[i][r.chain[i].size()-1]] - r.x[r.chain[i][0]]);
		r.distPBC(d);
		d = d.unit();
		dir += (dir*d > 0)? d : -d;
	}
	return dir.unit();
}

flt absPBC1(flt x, flt grid){
	while(x < 0)
		x += grid;
	while(x > grid)
		x -= grid;
	return x;
}

#define DISTR
#ifdef DISTR
flt countSmecticOrder(Restart const &r, flt periodRec) {
	flt periodStep = 0.05;
	flt shiftStep = 0.05;
	flt rcut = 3.f;
	int N = (int) round(r.box.sqr() / (M_PI * rcut * rcut));
	N = (N > 0)? N : 1;
	N *= 2;

	vect3 dir = countDirectorChains(r);
	vect3 planeX = vProduct(dir, vect3(1, 0, 0)).unit() * r.box.mod();
	vect3 planeY = vProduct(dir, planeX).unit() * r.box.mod();

	//find several distributions
	vector<vector<flt> > ncord;
	for(int n = 0; n < N; ++n) {
		vect3 r0 = r.box / 2 + planeX * ((flt) rand() / RAND_MAX - 1.f) + planeY * ((flt) rand() / RAND_MAX - 1.f);
		vector<flt> cord;
		for (int i = 0; i < r.nchains; ++i) {
			vect3 cm = centerOfMassOfChain(r, i);
			cm -= r0;
			r.distPBC(cm);
			flt c = dir * cm;
			if (sqrt(cm.sqr() - c * c) < rcut)
				cord.push_back(c);
		}
		if(cord.size() < 100) //not enough for statistics
			continue;
		ncord.push_back(cord);
	}

	// find period using first distribution
	vector<flt> cmin;
	flt stdevMin = -1;
	flt periodMin = -1;
	for (flt period = periodRec; period <= periodRec; period += periodStep)
		for (flt shift = 0; shift <= period; shift += shiftStep) {
			vector<flt> cp;
			for (int i = 0; i < ncord[0].size(); ++i)
				cp.push_back(absPBC1(ncord[0][i] + shift, period) / period);
			flt stdev, mean;
			meanAndStdev(cp, mean, stdev);
			if (stdev < stdevMin || stdevMin < 0) {
				cmin.swap(cp);
				stdevMin = stdev;
				periodMin = period;
			}
		}

	//add other distributions
	vector<flt> sumdist;
	for(int i = 0; i < cmin.size(); ++i)
		sumdist.push_back(cmin[i]);
	for(int i = 1; i < ncord.size(); ++i) {
		stdevMin = -1;
//		vector<flt> crand;
//		for(int j = 0; j < 5*ncord[i].size(); ++j)
//			crand.push_back((flt) rand() / (flt) RAND_MAX * periodMin);
		for (flt shift = 0; shift <= periodMin; shift += shiftStep) {
			vector<flt> cp;
			for (int j = 0; j < ncord[i-1].size(); ++j)
				cp.push_back(ncord[i-1][j]);
			for (int j = 0; j < ncord[i].size(); ++j)
				cp.push_back(absPBC1(ncord[i][j] + shift, periodMin) / periodMin);
//			for (int j = 0; j < crand.size(); ++j)
//				cp.push_back(absPBC1(crand[j] + shift, periodMin) / periodMin);
			flt stdev, mean;
			meanAndStdev(cp, mean, stdev);
			if (stdev < stdevMin || stdevMin < 0) {
				cmin = vector<flt>();
				for(int j = 0; j < ncord[i].size(); ++j)
					cmin.push_back(ncord[i][j] + shift);
				stdevMin = stdev;
			}
		}
		for(int j = 0; j < cmin.size(); ++j) {
			ncord[i][j] = cmin[j];
			sumdist.push_back(cmin[j]);
		}
	}

	//count dispersion
	stdevMin = -1;
	vector<flt> finDist;
	for(int i = 0; i < sumdist.size(); ++i)
		finDist.push_back(sumdist[i]);
	for(flt shift = 0.f; shift <= periodMin; shift += shiftStep) {
		vector<flt> cp;
		for (int i = 0; i < sumdist.size(); ++i)
			cp.push_back(absPBC1(sumdist[i] + shift, periodMin) / periodMin);
		flt stdev, mean;
		meanAndStdev(cp, mean, stdev);
		if (stdev < stdevMin || stdevMin < 0) {
			for(int i = 0; i < sumdist.size(); ++i)
				finDist[i] = cp[i];
			stdevMin = stdev;
		}
	}
	flt meanMean = 0;
	flt stdevMean = 0;
	meanAndStdev(finDist, meanMean, stdevMean);

//	cout << "period: " << periodMin << endl;
//	cout << "cords\tcords min\n";
//	for(int i = 0; i < cord.size(); ++i)
//		cout << cord[i] << '\t' << cmin[i] * periodMin << endl;
	stdevMean *= stdevMean;
	return 1.f - stdevMean * 12.f; // .289 - stdev of uniform distr
}
#else
flt countSmecticOrder(Restart const &r, flt periodRec) {
	flt periodStep = 0.05;
	flt shiftStep = 0.05;
	flt rcut = 3.f;
	int N = (int) round(r.box.sqr() / (M_PI * rcut * rcut));
	N = (N > 0)? N : 1;

	vect3 dir = countDirectorChains(r);
	vect3 planeX = vProduct(dir, vect3(1, 0, 0)).unit() * r.box.mod();
	vect3 planeY = vProduct(dir, planeX).unit() * r.box.mod();

	flt periodMean = 0;
	flt stdevMean = 0;
	int skipped = 0;
	for(int n = 0; n < N; ++n) {
		vect3 r0 = r.box / 2 + planeX * ((flt) rand() / RAND_MAX - 1.f) + planeY * ((flt) rand() / RAND_MAX - 1.f);
		vector<flt> cord;
		for (int i = 0; i < r.nchains; ++i) {
			vect3 cm = centerOfMassOfChain(r, i);
			cm -= r0;
			r.distPBC(cm);
			flt c = dir * cm;
			if (sqrt(cm.sqr() - c * c) < rcut)
				cord.push_back(c);
		}
		if(cord.size() < 100) { //not enough for statistics
			++skipped;
			continue;
		}

		vector<flt> cmin;
		flt stdevMin = -1;
		flt sdmReal;
		flt periodMin = -1;
		for (flt period = periodRec - .5f; period <= periodRec + .5f; period += periodStep) {
//			vector<flt> crand;
//			for(int i = 0; i < cord.size(); ++i)
//				crand.push_back((flt) rand() / (flt) RAND_MAX * period);
			for (flt shift = 0; shift <= period; shift += shiftStep) {
				vector<flt> cp;
				for (int i = 0; i < cord.size(); ++i) {
					cp.push_back(absPBC1(cord[i] + shift, period) / period);
//					cp.push_back(absPBC1(crand[i] + shift, period) / period);
				}
				flt stdev, mean;
				meanAndStdev(cp, mean, stdev);
				if (stdev < stdevMin || stdevMin < 0) {
					vector<flt> cpOrig;
					for (int i = 0; i < cord.size(); ++i)
						cpOrig.push_back(absPBC1(cord[i] + shift, period) / period);
					cmin.swap(cpOrig);
					meanAndStdev(cpOrig, mean, sdmReal);
					stdevMin = stdev;
					periodMin = period;
				}
			}
		}
		stdevMean += sdmReal;
		periodMean += periodMin;
	}
	N -= skipped;
	periodMean /= N;
	stdevMean /= N;

	cout << "period: " << periodMean << endl;
	cout << "skipped " << skipped << " of " << N+skipped << endl;
//	cout << "cords\tcords min\n";
//	for(int i = 0; i < cord.size(); ++i)
//		cout << cord[i] << '\t' << cmin[i] * periodMin << endl;
	stdevMean *= stdevMean;
	return 1.f - stdevMean * 12.f; // .289 - stdev of uniform distr
}
#endif

vect3 centerOfMassOfChain(Restart const &r, int i) {
	vect3 cm = r.x[r.chain[i][0]]; //center of mass of chain
	for (int j = 1; j < r.chain[i].size(); ++j) {
		vect3 dcm = r.x[r.chain[i][j]] - r.x[r.chain[i][0]];
		r.distPBC(dcm);
		cm += r.x[r.chain[i][0]] + dcm;
	}
	cm /= r.chain[i].size();
	return cm;
}

vector<vect3> countDisplacement(Restart &r1, Restart &r2, Units units, Space space, vect3 dir) {
	int n;
	if(units == ATOMS)
		n = r1.natoms;
	if(units == CHAINS)
		n = r1.nchains;
	vector<vect3> disp(n);
	
	for(int i = 0; i < n; ++i) {
		vect3 x1, x2;
		if(units == ATOMS) {
			x1 = r1.x[i];
			x2 = r2.x[i];
		}
		if(units == CHAINS) {
			x1 = centerOfMassOfChain(r1, i);
			x2 = centerOfMassOfChain(r2, i);
		}
		
		vect3 d = x2 - x1;
		r1.distPBC(d);
		if(space == LINE)
			d = (d*dir) * dir;
		if(space == PLANE)
			d = d - (d*dir) * dir;
		
		disp[i] = d;
	}
	return disp;
}







