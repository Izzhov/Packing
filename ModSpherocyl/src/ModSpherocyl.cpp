//============================================================================
// Name        : ModSpherocyl.cpp
// Author      : Kyle VanderWerf
// Version     :
// Copyright   : ORIGINAL CONTENT DO NOT STEAL
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include <ConNet.h>
#include <mm.h>
#include <Torus.h>
#include <Spherocyl.h>
#include <SpherocylDOF.h>
#include <HarmPot.h>
#include <HarmPotNbrList.h>
#include <SteepDesc.h>
#include <SteepDescAdapt.h>
#include <SecSteepDesc.h>
#include <FullSecSteepDesc.h>
#include <TimeRNGNormal.h>
#include <SingleMD.h>
#include <MultiMD.h>
#include <SimplexFltr.h>

#include "CompareHPs.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <boost/program_options.hpp>

#include <cmath>
#include <vector>

using namespace std;
namespace po = boost::program_options;

class mo{
public:
	template <class P, class B, class C>//potential,box,connet
	static void post_data(ofstream& dstream, const char * s,
				P p, B box, C cn){
		remove(s);
		dstream.open(s);
		//cout << s << endl;
		dstream << std::setprecision(10);
		dstream << "U: " << p.get_U() << "\n";
		dstream << "L: " << box.get_L() << "\n";
		dstream << "Packing Fraction: " << box.pack_frac() << "\n";
		dstream << "Seed: " << box.get_seed() << "\n";
		for(int i=0; i<box.get_N(); i++){
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_pos_coord(i,j) << " ";
			}
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_u_coord(i,j) << " ";
			}
			dstream << box.get_max_d(i) << " ";
			dstream << box.get_r(i) << " ";
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_1_pos_coord(i,j);
				if(j<(box.get_dim()-1)) dstream << " ";
			}
			dstream << endl;
		}
		cn.output_con(dstream);
		dstream.close();
	}
	template <class P, class B>//potential, box
	static void post_data(ofstream& dstream, const char * s,
				P p, B box){
		remove(s);
		dstream.open(s);
		//cout << s << endl;
		dstream << std::setprecision(10);
		dstream << "U: " << p.get_U() << "\n";
		dstream << "L: " << box.get_L() << "\n";
		dstream << "Packing Fraction: " << box.pack_frac() << "\n";
		dstream << "Seed: " << box.get_seed() << "\n";
		for(int i=0; i<box.get_N(); i++){
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_pos_coord(i,j) << " ";
			}
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_u_coord(i,j) << " ";
			}
			dstream << box.get_max_d(i) << " ";
			dstream << box.get_r(i) << " ";
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_1_pos_coord(i,j);
				if(j<(box.get_dim()-1)) dstream << " ";
			}
			dstream << endl;
		}
		dstream.close();
	}
	template <class T>
	static void post_box(ofstream& dstream, const char * s, Torus<T> box){
		remove(s);
		dstream.open(s);
		dstream << "L: " << box.get_L() << "\n";
		dstream << "Packing Fraction: " << box.pack_frac() << "\n";
		dstream << "Seed: " << box.get_seed() << "\n";
		for(int i=0; i<box.get_N(); i++){
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_pos_coord(i,j) << " ";
			}
			dstream << box.get_r(i) << " ";
			for(int j=0; j<box.get_dim(); j++){
				dstream << box.get_1_pos_coord(i,j);
				if(j<(box.get_dim()-1)) dstream << " ";
			}
			dstream << endl;
		}
		dstream.close();
	}
	static string to_string(int i){
		string result; ostringstream convert;
		convert << i; result = convert.str(); return result;
	}
	template<class P, class B>
	static void post_simp(int i, ofstream& dstream,
		P p, B box){
		string s = "pos_" + mo::to_string(i) + ".txt";
		const char * c = s.c_str();
		mo::post_data(dstream, c, p, box);
	}
};

int main(int argc, char **argv) {
	/*testing the distance algo
	//# particles, sizes, weights, dim, and init. packing frac
	int N=2; int dim=2; double L=6.487363085;
	vector<vector<gsl_vector*> > szs;
	vector<vector<gsl_vector*> > locs;
	gsl_vector * r1 = gsl_vector_alloc(1); gsl_vector_set(r1,0,1);
	gsl_vector * r2 = gsl_vector_alloc(1); gsl_vector_set(r2,0,1);
	gsl_vector * a1 = gsl_vector_alloc(1); gsl_vector_set(a1,0,1.2);
	gsl_vector * a2 = gsl_vector_alloc(1); gsl_vector_set(a2,0,1.2);
	vector<gsl_vector*> s1; s1.push_back(r1); s1.push_back(a1);
	vector<gsl_vector*> s2; s2.push_back(r2); s2.push_back(a2);
	szs.push_back(s1); szs.push_back(s2);
	gsl_vector * l1 = gsl_vector_alloc(dim);
	gsl_vector_set(l1,0,0.3430112694); gsl_vector_set(l1,1,0.3047892792);
	gsl_vector * u1 = gsl_vector_alloc(dim);
	gsl_vector_set(u1,0,0.6358887663); gsl_vector_set(u1,1,0.7717807181);
	gsl_vector * l2 = gsl_vector_alloc(dim);
	gsl_vector_set(l2,0,0.08472802131); gsl_vector_set(l2,1,0.4786794158);
	gsl_vector * u2 = gsl_vector_alloc(dim);
	gsl_vector_set(u2,0,0.6756162572); gsl_vector_set(u2,1,0.737253466);
	vector<gsl_vector*> loc1; loc1.push_back(l1); loc1.push_back(u1);
	vector<gsl_vector*> loc2; loc2.push_back(l2); loc2.push_back(u2);
	locs.push_back(loc1); locs.push_back(loc2);
	Torus<SpherocylDOF> box(N,L,locs,szs,dim);
	cout << setprecision(12);
	cout << std::sqrt(box.ell2(0,1,1,0)) << endl;
	cout << std::sqrt(box.ell2(0,1,1,1)) << endl;
	ConNet<Torus<SpherocylDOF> > cn(box);
	cn.print_con();
	*/

	po::options_description desc("Allowed options");
	desc.add_options()
			("seed", po::value<int>(), "set seed")
			("filename", po::value<string>(), "set filename")
			("ar", po::value<double>(), "set aspect ratio")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);

	//# particles, sizes, weights, dim, and init. packing frac
	int N=6; int dim=2; double phi=0.01;
	vector<vector<gsl_vector*> > szs;
	vector<double> wts;
	vector<gsl_vector*> sz1;
	vector<gsl_vector*> sz2;
	gsl_vector * r1 = gsl_vector_alloc(1);
	gsl_vector * a1 = gsl_vector_alloc(1);
	gsl_vector * r2 = gsl_vector_alloc(1);
	gsl_vector * a2 = gsl_vector_alloc(1);
	gsl_vector_set(r1,0,1.0); gsl_vector_set(a1,0,1.2);
	gsl_vector_set(r2,0,1.4); gsl_vector_set(a2,0,1.2);
	if(vm.count("ar")){
		gsl_vector_set(a1,0,vm["ar"].as<double>());
		gsl_vector_set(a2,0,vm["ar"].as<double>());
	}
	sz1.push_back(r1); sz1.push_back(a1);
	sz2.push_back(r2); sz2.push_back(a2);
	szs.push_back(sz1); szs.push_back(sz2);
	wts.push_back(1.0); wts.push_back(1.0);
	Torus<SpherocylDOF> box(N,szs,wts,dim,phi);
	//Torus<SpherocylDOF> box2(N,szs,wts,dim,phi,1435762516);
	//1433971798 000@3.3 (wha happun?)
	//1433972169,1433972110 (some 3.3's)
	//1435762516 new bidi1 (oldbidi3)
	//1435691590 doesn't seem stable (turns out I was missing one)
	if(vm.count("seed")){
		Torus<SpherocylDOF> box2(N,szs,wts,dim,phi,vm["seed"].as<int>());
		box = box2;
	}
	HarmPot<Torus<SpherocylDOF> > pot(box);
	//HarmPotNbrList<Torus<SpherocylDOF>,SpherocylDOF > pot2(box2);
	FullSecSteepDesc<HarmPot<Torus<SpherocylDOF> >, Torus<SpherocylDOF> > min(pot, box);
	//SecSteepDesc<HarmPotNbrList<Torus<SpherocylDOF>,SpherocylDOF >, Torus<SpherocylDOF> > min2(pot2, box2);
	//CompareHPs min(pot,box,pot2,box2);

	ofstream dstream;

	bool fltt = min.exp_minimize(7);
	//bool fltt2 = min2.minimize();
	//bool fltt = false;

	cout << box.get_seed() << endl;

	//cout << "PackFrac: " << box.pack_frac() << endl;

	ConNet<Torus<SpherocylDOF> > cn(box);
	//ConNet<Torus<SpherocylDOF> > cn2(box2);
	if(vm.count("filename")){
		string name = vm["filename"].as<string>();
		const char * charname = name.c_str();
		mo::post_data(dstream, charname, pot, box,cn);
	}
	else mo::post_data(dstream, "fin_pos.txt", pot, box,cn);
	//mo::post_data(dstream, "fin_pos2.txt", pot2, box2,cn2);


	cout << "Final Pressure: " << pot.pressure() << endl;
	dstream.open("aspect_data.txt", ofstream::out | ofstream::app);
	dstream << box.get_seed() << " " << box.get_shape(0).get_a();
	dstream << " " << cn.num_contacts() << " " << fltt << " ";
	dstream << min.eps_reached() << endl;
	dstream.close();

	//cn.print_con();

	/*testing simplexfltr
	int d=3; int sym=2;
	vector<gsl_vector*> fis; vector<gsl_vector*> flocs;
	gsl_vector * v1 = gsl_vector_alloc(2);
	gsl_vector_set(v1,0,1); gsl_vector_set(v1,1,-1);
	gsl_vector * v2 = gsl_vector_alloc(2);
	gsl_vector_set(v2,0,-3); gsl_vector_set(v2,1,-1);
	gsl_vector * v3 = gsl_vector_alloc(2);
	gsl_vector_set(v3,0,0); gsl_vector_set(v3,1,1);
	gsl_vector * v4 = gsl_vector_alloc(2);
	gsl_vector_set(v4,0,0); gsl_vector_set(v4,1,1);
	fis.push_back(v1); fis.push_back(v2);
	fis.push_back(v3); fis.push_back(v4);
	gsl_vector * w1 = gsl_vector_alloc(2);
	gsl_vector_set(w1,0,-2); gsl_vector_set(w1,1,0);
	gsl_vector * w2 = gsl_vector_alloc(2);
	gsl_vector_set(w2,0,2); gsl_vector_set(w2,1,0);
	gsl_vector * w3 = gsl_vector_alloc(2);
	gsl_vector_set(w3,0,2); gsl_vector_set(w3,1,0);
	gsl_vector * w4 = gsl_vector_alloc(2);
	gsl_vector_set(w4,0,-2); gsl_vector_set(w4,1,0);
	flocs.push_back(w1); flocs.push_back(w2);
	flocs.push_back(w3); flocs.push_back(w4);
	SimplexFltr flt_test(fis,flocs,d,sym);
	bool isit = flt_test.is_floater();
	std::cout << isit << std::endl;
	*/
	return 0;
}
