//Algorithm.cpp

#include "Algorithm.h"

Algorithm::Algorithm(){}

Algorithm::~Algorithm(){
    for(unsigned int i=0; i<logs.size(); i++){
      delete logs[i];
    }
}

