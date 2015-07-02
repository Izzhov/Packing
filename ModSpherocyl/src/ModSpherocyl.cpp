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
#include <HarmPot.h>
#include <SteepDesc.h>
#include <TimeRNGNormal.h>
#include <SingleMD.h>
#include <MultiMD.h>

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
	int N=2; int dim=2; double L=10.0;
	vector<vector<gsl_vector*> > szs;
	vector<vector<gsl_vector*> > locs;
	gsl_vector * r1 = gsl_vector_alloc(1); gsl_vector_set(r1,0,1);
	gsl_vector * r2 = gsl_vector_alloc(1); gsl_vector_set(r2,0,1);
	gsl_vector * a1 = gsl_vector_alloc(1); gsl_vector_set(a1,0,2);
	gsl_vector * a2 = gsl_vector_alloc(1); gsl_vector_set(a2,0,2);
	vector<gsl_vector*> s1; s1.push_back(r1); s1.push_back(a1);
	vector<gsl_vector*> s2; s2.push_back(r2); s2.push_back(a2);
	szs.push_back(s1); szs.push_back(s2);
	gsl_vector * l1 = gsl_vector_alloc(dim);
	gsl_vector_set(l1,0,0.5); gsl_vector_set(l1,1,0.5);
	gsl_vector * u1 = gsl_vector_alloc(dim);
	gsl_vector_set(u1,0,1); gsl_vector_set(u1,1,0);
	gsl_vector * l2 = gsl_vector_alloc(dim);
	gsl_vector_set(l2,0,0.7); gsl_vector_set(l2,1,0.6);
	gsl_vector * u2 = gsl_vector_alloc(dim);
	gsl_vector_set(u2,0,0); gsl_vector_set(u2,1,1);
	vector<gsl_vector*> loc1; loc1.push_back(l1); loc1.push_back(u1);
	vector<gsl_vector*> loc2; loc2.push_back(l2); loc2.push_back(u2);
	locs.push_back(loc1); locs.push_back(loc2);
	Torus<Spherocyl> box(N,L,locs,szs,dim);
	HarmPot<Torus<Spherocyl> > pot(box);
	SteepDesc<HarmPot<Torus<Spherocyl> >, Torus<Spherocyl> > min(pot, box);
	min.next();
	cout << box.get_1_pos_coord(0,0) << " " << box.get_1_pos_coord(0,1) << endl;
	cout << box.get_u_coord(0,0) << " " << box.get_u_coord(0,1) << endl;
	cout << box.get_1_pos_coord(1,0) << " " << box.get_1_pos_coord(1,1) << endl;
	cout << box.get_u_coord(1,0) << " " << box.get_u_coord(1,1) << endl;
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
	Torus<Spherocyl> box(N,szs,wts,dim,phi);
	//1435762516 new bidi1 (oldbidi3)
	//1435691590 doesn't seem stable (turns out I was missing one)
	if(vm.count("seed")){
		Torus<Spherocyl> box2(N,szs,wts,dim,phi,vm["seed"].as<int>());
		box = box2;
	}
	HarmPot<Torus<Spherocyl> > pot(box);
	SteepDesc<HarmPot<Torus<Spherocyl> >, Torus<Spherocyl> > min(pot, box);

	ofstream dstream;

	bool fltt = min.exp_minimize(10);

	//cout << box.get_seed() << endl;

	//cout << "PackFrac: " << box.pack_frac() << endl;

	ConNet<Torus<Spherocyl> > cn(box);
	if(vm.count("filename")){
		string name = vm["filename"].as<string>();
		const char * charname = name.c_str();
		mo::post_data(dstream, charname, pot, box,cn);
	}
	/*
	else mo::post_data(dstream, "fin_pos.txt", pot, box,cn);

	cout << "Final Pressure: " << pot.pressure() << endl;
	cout << box.ell2(0,2,0) << endl;
	*/
	dstream.open("aspect_data.txt", ofstream::out | ofstream::app);
	dstream << box.get_seed() << " " << box.get_shape(0).get_a();
	dstream << " " << cn.num_contacts() << " " << fltt << endl;
	dstream.close();
	return 0;
}
