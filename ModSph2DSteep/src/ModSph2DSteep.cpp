//============================================================================
// Name        : ModSph2DSteep.cpp
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

#include <mm.h>
#include <Torus.h>
#include <Sphere.h>
#include <HarmPot.h>
#include <SteepDesc.h>
#include <SphereConNet.h>
#include <SphDynMatHarmOsc.h>
#include <TimeRNGNormal.h>
#include <SingleMD.h>

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
	template <class P, class B>//potential, box
	static void post_data(ofstream& dstream, const char * s,
				P p, B box, SphereConNet cn){
		remove(s);
		dstream.open(s);
		dstream << std::setprecision(10);
		dstream << "U: " << p.get_U() << "\n";
		dstream << "L: " << box.get_L() << "\n";
		dstream << "Packing Fraction: " << box.pack_frac() << "\n";
		dstream << "Seed: " << box.get_seed() << "\n";
		dstream << "Number of Contacts: " << cn.num_contacts() << "\n";
		dstream << "Desired # of Contacts: " << cn.desired_contacts() << "\n";
		std::vector<int> whfl = cn.which_floaters();
		dstream << "Which Floaters: ";
		for(unsigned int i=0;i<whfl.size();i++) dstream << whfl.at(i) << " ";
		dstream << "\n";
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
		cn.output_con(dstream);
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
		P p, B box, SphereConNet cn){
		string s = "pos_" + mo::to_string(i) + ".txt";
		const char * c = s.c_str();
		mo::post_data(dstream, c, p, box, cn);
	}
};

int main(int argc, char **argv) {

	po::options_description desc("Allowed options");
	desc.add_options()
			("seed", po::value<int>(), "set seed")
			("filename", po::value<string>(), "set filename")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);

	//# particles, sizes, weights, dim, and init. packing frac
	int N = 6;
	vector<vector<gsl_vector*> > szs;
	vector<double> wts;
	int dim = 2;
	double phi = 0.01;
	vector<gsl_vector*> sz1; vector<gsl_vector*> sz2;
	gsl_vector * size1 = gsl_vector_alloc(1);
	gsl_vector * size2 = gsl_vector_alloc(1);
	gsl_vector_set(size1, 0, 1.0); gsl_vector_set(size2, 0, 1.4);
	sz1.push_back(size1); sz2.push_back(size2);
	szs.push_back(sz1); szs.push_back(sz2);
	wts.push_back(1.0); wts.push_back(1.0);
	Torus<Sphere> box(N,szs,wts,dim,phi);
	//1434000791 for 6 is null and PF not 8
	//1434653545
	//1435088168 for floater (needs 10^-12 to be accurate)
	if(vm.count("seed")){
		Torus<Sphere> box2(N,szs,wts,dim,phi,vm["seed"].as<int>());
		box = box2;
	}
	HarmPot<Torus<Sphere> > pot(box);
	SteepDesc<HarmPot<Torus<Sphere> >, Torus<Sphere> > min(pot, box);

	ofstream dstream;

	min.minimize(5);

	cout << box.get_seed() << endl;

	SphereConNet cn(box,true);//will auto remove floaters

	if(vm.count("filename")){
		string name = vm["filename"].as<string>();
		const char * charname = name.c_str();
		mo::post_data(dstream, charname, pot, box, cn);
	}
	else mo::post_data(dstream, "fin_pos.txt", pot, box, cn);
	/*
	cout << "Number of Contacts: " << cn.num_contacts() << endl;
	cout << "Desired Contacts: " << cn.desired_contacts() << endl;
	cout << "Final Pressure: " << pot.pressure() << endl;

	SphDynMatHarmOsc sdmho(box,pot,cn);

	remove("dynmat.txt");
	dstream.open("dynmat.txt");
	dstream << setprecision(12);
	for(int i=0; i<12; i++){
		for(int j=0; j<12; j++){
			dstream << gsl_matrix_get(sdmho.get_mat(),i,j);
			if(j<12) dstream << " ";
		}
		dstream << endl;
	}

	TimeRNGNormal trng;
	SingleMD<HarmPot<Torus<Sphere> >, Torus<Sphere> > smd(box,0.1,trng);

	cout << "all finished." << endl;

	cn.print_con();
	*/
	return 0;
}
