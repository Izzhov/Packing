//============================================================================
// Name        : Sphere2DSteepest.cpp
// Author      : Kyle VanderWerf
// Version     :
// Copyright   : ORIGINAL CONTENT DO NOT STEAL
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>

#include <Sphere2D.h>
#include <SphereSizes.h>
#include <TimeRNG01.h>
#include <Torus2D.h>
#include <mymath.h>
#include <Potential2DNoNbr.h>
#include <SphereConNet2D.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <boost/program_options.hpp>

#include <stdio.h>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>

#include <gsl/gsl_linalg.h>

using namespace std;
namespace po = boost::program_options;

class mm{
public:
	static void post_data(ofstream& dstream, const char * s,
			Potential2DNoNbr p, Torus2D box, SphereConNet2D cn){
		remove(s);
		dstream.open(s);
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
			dstream << box.get_x(i) << " " << box.get_y(i);
			dstream << " " << box.get_r(i) << " ";
			dstream << box.get_1x(i) << " " << box.get_1y(i) << "\n";
		}
		dstream.close();
	}
	static string to_string(int i){
		string result; ostringstream convert;
		convert << i; result = convert.str(); return result;
	}
	static void post_simp(int i, ofstream& dstream,
			Potential2DNoNbr p, Torus2D box, SphereConNet2D cn){
		string s = "pos_" + mm::to_string(i) + ".txt";
		const char * c = s.c_str();
		mm::post_data(dstream, c, p, box, cn);
	}
};

int main(int argc, char **argv) {
	//system("rm *.txt");
	/*
	 * set: # spheres
	 *      packing fraction for initial area A
	 *      sphere sizes
	 */
	// random seed

	po::options_description desc("Allowed options");
	desc.add_options()
			("seed", po::value<int>(), "set seed")
			("filename", po::value<string>(), "set filename")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);

	int N = 16;
	double phi = 0.01;
	std::vector<double> sizes; std::vector<double> wts;
	sizes.push_back(1.0); sizes.push_back(1.4);
	wts.push_back(1.0); wts.push_back(1.0);
	Torus2D box(N, sizes, wts);
	if(vm.count("seed")){
		Torus2D box2(N,sizes,wts,vm["seed"].as<int>());
		box = box2;
	}
	box.find_set_L(phi);

	Potential2DNoNbr p(box);
	ofstream dstream;

	//mm::post_simp(0,dstream,p,box);


	//check when they touch
	//int count = 0;
	for(int i=1; i<100000*10*100; i++){
		p.full_next();
		/* for when they touch
		if (p.get_U() - p.get_P()*box.get_L()*box.get_L()>0.000001){
			cout << "i=" << i << "; " << "U-LU=";
			cout << p.get_U()-p.get_P()*box.get_L()*box.get_L() << endl;
			mm::post_simp(i,dstream,p,box);
			count++; if(count>100) break;
		}
		*/
		if(i%1000==0) if(p.is_done()) break;
		//if(i%(100000)==0) mm::post_simp(i,dstream,p,box);
	}

	/*epscheck
	double e;
	for(int i = 1; i<100000000; i++){
		e=p.next_epscheck();
		p.next_rest();
		if(e>0.00011 || e<0.00009){
			cout << "i="; cout << i; cout << "; e=";
			cout << e << endl;
			mm::post_simp(i,dstream,p,box);
		}
		else if (i % 10000000 == 0) mm::post_simp(i, dstream, p, box);
	}
	*/

	/*
	for(int j = 2; j<10; j++){
		for(int i=1; i<100000; i++) p.full_next();
		s = "pos_" + mm::to_string(j) + ".txt";
		c = s.c_str();
		mm::post_data(dstream, c, p, box);
	}
	*/
	//check with pg. 5-27-15
	/*
	int N = 2; double L = 10.0;
	vector<vector<double> > locs;
	vector<double> a; a.push_back(0.6); a.push_back(0);
	vector<double> b; b.push_back(0.5); b.push_back(0);
	locs.push_back(a); locs.push_back(b);
	Torus2D box(N,L,locs);

	Potential2DNoNbr p(box);
	ofstream dstream;

	mm::post_simp(0,dstream,p,box);

	cout << p.dxi_ell(0,1,2) << endl;

	p.full_next();
	*/

	SphereConNet2D cn(box);//will auto remove floaters

	if(vm.count("filename")){
		string name = vm["filename"].as<string>();
		const char * charname = name.c_str();
		mm::post_data(dstream, charname, p, box, cn);
	}
	else mm::post_data(dstream, "fin_pos.txt", p, box,cn);


	//cn.print_con();
	//cout << box.ell(1,8,1) << endl;
	//cout << box.ell(1,13,2) << endl;

	//TimeRNG01 trng;
	//cout << trng.num() << endl;

	return 0;
}
