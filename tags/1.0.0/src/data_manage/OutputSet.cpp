#include "OutputSet.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
using namespace std;

namespace data_manage
{

///
/// outputs designated set
///
void OutputSet::outputSet(string name, Dataset& set){

    ofstream outfile;
    outfile.open(name.c_str(), ios::out);    
    outfile << set;
    outfile.close();
   
}

///
/// outputs set and adds cv number to name
///
void OutputSet::outputCV(string name, Dataset& set, int cv){

  stringstream ss;
  ss << name << "." << cv << ".txt";
  outputSet(ss.str(), set);

}

}
