//Algorithm.cpp

#include "Algorithm.h"

Algorithm::Algorithm(){
    fitness_name=" ";
}

Algorithm::~Algorithm(){
    for(unsigned int i=0; i<logs.size(); i++){
      delete logs[i];
    }
}

