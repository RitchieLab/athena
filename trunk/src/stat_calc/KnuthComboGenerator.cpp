#include "KnuthComboGenerator.h"
using namespace std;

namespace stat{

///
/// Constructor initializes parameters
///
KnuthComboGenerator::KnuthComboGenerator(){
  initialize();
}

///
/// Destructor
///
KnuthComboGenerator::~KnuthComboGenerator(){
// 	if(c != NULL)
// 		delete [] c;
}

///
/// Initializes variables for generator
/// @return none
///
void KnuthComboGenerator::initialize(){
//   c = NULL;
  ComboInterval = 1000;
  AlreadyStarted = false;
  NumLoci = 0;
  ComboStart=0;
  ComboEnd=0;
}


///
/// Creates initial state 
///
void KnuthComboGenerator::initialize_state(){
  kdec = ComboEnd;
  x = 0;
  j = kdec;

// 	if(c != NULL)
// 		delete c;
//   c = new int[ComboEnd+3];
	c.resize(ComboEnd+3);

  for(int i=1; i <= kdec; i++){
    c[i] = i;
  }
  c[kdec+1] = NumLoci+1;
  c[kdec+2] = 0;
}


///
/// Generates combinations in amount passed
/// @param new_interval Combinations to create
/// 
bool KnuthComboGenerator::GenerateCombinations(int new_interval){
  int old_interval = getComboInterval();
  SetComboInterval(new_interval);
  bool done = GenerateCombinations();
  SetComboInterval(old_interval);
  return done;
}

///
/// Uses modified code by Donald Knuth and Glenn C. Rhoads to generate
/// all possible combinations from COMBOSTART to COMBOEND in intervals specified in COMBOINTERVAL
/// That is, N = COMBOEND, C = COMBOSTART. Summation of N choose C  -->  N choose N.
/// Code is confusing but fast
/// @return true if all combinations completed
///
bool KnuthComboGenerator::GenerateCombinations(){
  // Clear the ComboList!!! Very important, or Memory will overflow
  ComboList.clear();

  counter = 0;  // counter for the number of combinations created

  if(AlreadyStarted)   // If combination generator has already been started
  {
    counter--;
    goto resume;    // Ugh.. a GOTO... but it works!
  }
  else{
    kdec = ComboEnd;
    x = 0;
    j = kdec;
    
    if(NumLoci == ComboEnd){
      j=0;
    }
    
//     c = new int[ComboEnd+3];
    c.resize(ComboEnd+3);

    for(int i=1; i <= kdec; i++)
    {
      c[i] = i;
    }
    c[kdec+1] = NumLoci+1;
    c[kdec+2] = 0;
  }

  AlreadyStarted = true;  // If it wasn't already started, it is now

  init:

      if(kdec == (ComboStart-1))
      {
//         delete [] c;
        c.clear();
        AlreadyStarted = false;
        counter--;
        return true;
      }

      for (int i=1; i <= kdec; i++)
      {  c[i] = i;  }
      c[kdec+1] = NumLoci+1;
      c[kdec+2] = 0;
      j = kdec;
      
      if(NumLoci == ComboEnd){
        j=0;
      }
      
      // Add a new combination to the ComboList
      ComboList.push_back(std::vector <unsigned int>());

  visit:
      for (int i=kdec; i >= 1; i--)
      {
    // Add an element to the new combination
    ComboList[counter].push_back(c[i]);
      }


      if(counter >= ComboInterval)  // If you exceed the interval limit
      {
          return false;
      }

  resume:
            // Add a new combination to the ComboList
      ComboList.push_back(std::vector <unsigned int>());

      counter++;
      if (j > 0) {x = j+1; goto incr;}

      if (c[1] + 1 < c[2])
         {
         c[1] += 1;
         goto visit;
         }

      j = 2;

   do_more:
      c[j-1] = j-1;
      x = c[j] + 1;
      if (x == c[j+1]) {j++; goto do_more;}

      if (j > kdec)
      {
          kdec--;

    // Remove the last empty combination before quitting
          ComboList.pop_back();
          goto init;
      }

   incr:
      c[j] = x;
      j--;
      goto visit;
}

///
/// returns total number of combinations in set
/// @return number of combinations to construct
///
double KnuthComboGenerator::calc_combinations(){
  double total = 0.0;
  
  unsigned int endFact, n = NumLoci;
  double k, numerator, k_fact, result;
  unsigned int i;
  
  for(int combsize=ComboStart; combsize <= ComboEnd;
    combsize++){
    k = combsize;
    numerator = double(n);
    endFact = n-combsize;
    // calculate numerator
    for(i=NumLoci-1; i>endFact; i--){
      numerator *= i;
    }

    k_fact = combsize;
    for(i=combsize-1; i >=1; i--){
      k_fact *= i;
    }

    result = numerator / k_fact;
    total += result;
  }
  return total;
}

///
///Sets the interval for generating combinations.  Used so that
///don't excessively use memory 
//@param cmbInterval interval to use in producing combinations
//@return none
///
void KnuthComboGenerator::SetComboInterval(int cmbInterval){
  ComboInterval = cmbInterval;
}


///
///Sets the start combination number and ending combination
///@param  combStart minimum size of combination
///@param  combEnd maximum size of combination
///@return  none
///
void KnuthComboGenerator::ComboEnds(int combStart, int combEnd){
  ComboStart = combStart;
  ComboEnd = combEnd;
//  initialize();
}

///
///Sets total number of loci for generator
///@param nLoci number of loci
///@return none
///
void KnuthComboGenerator::SetLoci(int nLoci){
  NumLoci = nLoci;
}

}

