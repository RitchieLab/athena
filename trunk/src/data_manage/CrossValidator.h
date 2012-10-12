/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CROSSVALIDATOR_H_
#define CROSSVALIDATOR_H_

#include <cstdlib>

#include "CVSet.h"
#include "Dataholder.h"
#include "random_func.h"

namespace data_manage
{


class TRandom
{
  public:
  int operator()(int n)
  {
//      return rand() % n;
	return rand_uint() % n;
  }

};

///
///  Performs crossvalidation on a dataholder to return datasets with the
///

class CrossValidator
{
public:
  CrossValidator();
  ~CrossValidator();

  /// splits data into testing and training sets for the indicated number of intervals
  CVSet split_data(unsigned int num_crossval, Dataholder* holder);
  
  /// save splits with individual IDs
  void save_splits(std::string filename);
  
  /// load splits with individual IDs
  CVSet load_splits(std::string filename, Dataholder* holder);
  
//   friend std::ostream& operator<<( std::ostream& os, const CrossValidator& cv );
//   friend std::istream& operator>>( std::istream& is, CrossValidator& cv );

private:

  void distribute_inds(unsigned int num_splits, std::vector<Individual*>& inds,
      std::vector<std::vector<Individual*> >& splits);

  void shuffle_inds(std::vector<Individual*> & inds);

  void status_bin(Dataholder* holder, std::vector<Individual*>& affected,
      std::vector<Individual*>& unaffected);
      
  CVSet split_by_num(Dataholder* holder);    
  
  CVSet create_set(unsigned int num_cv, Dataholder* holder);

  vector<vector<Individual*> > splits;
  
  Dataholder* ref_holder;

};

}

#endif /*CROSSVALIDATOR_H_*/
