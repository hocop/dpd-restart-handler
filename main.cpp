#include "Analyser.h"
#include <dirent.h>
#include <sstream>
#include <iomanip>

bool directoryExists( const char* pzPath );

int main() {
//	Restart r = importGav("/home/rusl/Data/lc/comp_RC8eqrod/r7c2-0.40/restart20.dat");
//	exportNormal(r, "/home/rusl/Data/lc/comp_RC8eqrod/r7c2-0.40/restart20.ndat");
//	Restart r2 = importNormal("/home/rusl/Data/lc/comp_RC8eqrod/r7c2-0.40/restart20.ndat");
//	exportSmall(r2, "/home/rusl/Data/lc/comp_RC8eqrod/r7c2-0.40/restart20.sdat");
//	Restart r3 = importSmall("/home/rusl/Data/lc/comp_RC8eqrod/r7c2-0.40/restart20.sdat");
//	exportGav(r3, "/home/rusl/Data/lc/comp_RC8eqrod/r7c2-0.40/restart20-test.dat");
//	return 0;
	
//	string directory = "/home/rusl/Data/lc/comp/r7c_den3/";
//	string directory = "/media/rusl/MyPassport/comp/r7c_den4/";
	string directory = "/run/user/1000/gvfs/sftp:host=lom2/home/vurdizm/_scratch/students/bainazarov/lc_phases/comp_s/r7c_den4/";
	
	vector<string> name =      parseLine("r7c0 r7c1 r7c2 r7c3 r7c4 r7c5 r7c6 r7c7");
	vector<string> periodApr = parseLine("3.5  4.0  4.5  5.0  5.5  6.0  6.5  7.0 ");
//	vector<string> periodApr = parseLine("4.667  5.333  6.0  6.667  7.333  8.0  8.667  9.333 ");
//	vector<string> file = parseLine("restart10.dat restart11.dat restart12.dat restart13.dat restart14.dat restart15.dat restart16.dat restart17.dat restart18.dat restart19.dat restart20.dat restart21.dat");
	vector<string> file = parseLine("restart10.sdat restart11.sdat restart12.sdat restart13.sdat restart14.sdat restart15.sdat restart16.sdat restart17.sdat restart18.sdat restart19.sdat restart20.sdat restart21.sdat");
//	vector<string> file = parseLine("restart20.dat");
	
	vector<string> temp;
	for(flt t1 = 0.30, t2 = 0.305; t1 <= 1.00; t1 += 0.01, t2 += 0.01) {
		stringstream stream1, stream2;
		stream1 << fixed << setprecision(2) << t1;
		temp.push_back(stream1.str());
		stream2 << fixed << setprecision(3) << t2;
		temp.push_back(stream2.str());
	}
	
	int rmin = 3;
	int rmax = 20;
	
	vector<int> type;
	type.push_back(3);
	
	vector<vector<flt> > nematic;
	vector<flt> smectic;
	vector<bool> tempFound;
	for(int i = 0; i < name.size(); ++i) {
		flt period = (flt) atof(periodApr[i].c_str());
		for (int j = 0; j < temp.size(); ++j) {
			string dirtemp = directory + name[i] + '-' + temp[j];
			bool tfound = directoryExists(dirtemp.c_str());
			tempFound.push_back(tfound);
			
			if (tfound) {
				vector<vector<flt> > nem;
				vector<flt> sm;
				for (int k = 0; k < file.size(); ++k) {
					string s = dirtemp + '/' + file[k];
					cout << "analysing " << s << endl;
					Restart r0 = importSmall(s);
					Restart r = selectByType(r0, type);
					nem.push_back(countNematicOrderChains(r, rmin, rmax));
					sm.push_back(countSmecticOrder(r, period));
				}
				// averaging smectic
				flt smMean = 0;
				for (int k = 0; k < file.size(); ++k)
					smMean += sm[k];
				smMean /= file.size();
				smectic.push_back(smMean);
				
				//averaging nematic
				vector<flt> nemMean(nem[0].size());
				for (int k = 0; k < nemMean.size(); ++k)
					nemMean[k] = 0.f;
				for (int k = 0; k < file.size(); ++k)
					for (int m = 0; m < nemMean.size(); ++m)
						nemMean[m] += nem[k][m];
				for (int k = 0; k < nemMean.size(); ++k)
					nemMean[k] /= file.size();
				
				nematic.push_back(nemMean);
			} else {
				smectic.push_back(-1.f);
				vector<flt> nem(rmax + 1 - rmin);
				for (int k = 0; k < nem.size(); ++k)
					nem[k] = -1.f;
				nematic.push_back(nem);
			}
		}
	}
	
	//writing output
	// nematic
//	cout << "Nematic order S(R)\n";
//	for(int i = 0; i < name.size(); ++i) {
////		 names of columns
//		cout << name[i];
//		for (int j = 0; j < temp.size(); ++j)
//			cout << '\t' << temp[j];
//
////		 values
//		for(int k = 0; k < nematic[0].size(); ++k) {
//			cout << '\n';
//			cout << rmin + k;
//			for (int j = 0; j < temp.size(); ++j)
//				if(tempFound[i*temp.size() + j])
//					cout << '\t' << nematic[i*temp.size() + j][k];
//				else
//					cout << '\t' << '~';
//		}
//		cout << "\n\n";
//	}
	//nematic total
	cout << "Nematic order at max R\n";
	cout << "T";
	for(int i = 0; i < name.size(); ++i)
		cout << "\t" << name[i];
	cout << endl;
	for(int i = 0; i < temp.size(); ++i) {
		bool prntThis = false;
		for(int j = 0; j < name.size(); ++j)
			if(tempFound[j*temp.size() + i]) {
				prntThis = true;
				break;
			}
		if(prntThis) {
			cout << temp[i];
			for (int j = 0; j < name.size(); ++j)
				if (tempFound[j * temp.size() + i])
					cout << '\t' << nematic[j * temp.size() + i][rmax - rmin - 1];
				else
					cout << '\t' << '~';
			cout << endl;
		}
	}
	//smectic
	cout << "\n\nSmectic order\n";
	cout << "T";
	for(int i = 0; i < name.size(); ++i)
		cout << "\t" << name[i];
	cout << endl;
	for(int i = 0; i < temp.size(); ++i) {
		bool prntThis = false;
		for(int j = 0; j < name.size(); ++j)
			if(tempFound[j*temp.size() + i]) {
				prntThis = true;
				break;
			}
		if(prntThis) {
			cout << temp[i];
			for (int j = 0; j < name.size(); ++j)
				if (tempFound[j * temp.size() + i])
					cout << '\t' << smectic[j * temp.size() + i];
				else
					cout << '\t' << '~';
			cout << endl;
		}
	}
	
	//smectic and nematic
	for(int i = 0; i < name.size(); ++i) {
		cout << "\n\nSmectic and nematic order for " << name[i] << endl;
		cout << name[i] << "\tnematic" << "\tsmectic\n";
		for(int j = 0; j < temp.size(); ++j)
			if(tempFound[i*temp.size() + j])
				cout << temp[j] << "\t" << nematic[i*temp.size() + j][rmax-rmin] << '\t' << smectic[i*temp.size() + j] << endl;
		cout << "\n\n";
	}
	
	return 0;
}




bool directoryExists( const char* pzPath )
{
	if ( pzPath == NULL) return false;
	
	DIR *pDir;
	bool bExists = false;
	
	pDir = opendir (pzPath);
	
	if (pDir != NULL)
	{
		bExists = true;
		(void) closedir (pDir);
	}
	
	return bExists;
}









