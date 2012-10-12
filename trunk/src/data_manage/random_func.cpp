#include "random_func.h"
#include "../random/mersenne_twister.h"
#include <iostream>
#include <fstream>
using namespace std;

namespace data_manage{

static Random::MTRand mersenne_twister;

void rand_seed(unsigned int seed){
	mersenne_twister.seed(seed);
}

float rand_float(){
	return mersenne_twister.rand();
}

unsigned int rand_uint(){
	return mersenne_twister.randInt();
}


void rand_load_state(string filename){
	ifstream is;
	is.open(filename.c_str());
	is >> mersenne_twister;
	is.close();
}

void rand_save_state(string filename){
	ofstream os(filename.c_str(), ios::out);
	os << mersenne_twister;
	os.close();
}

}