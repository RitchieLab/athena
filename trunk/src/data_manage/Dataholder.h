#ifndef DATAHOLDER_H_
#define DATAHOLDER_H_

#include "Individual.h"
#include "DataExcept.h"
#include <string>

using namespace std;

namespace data_manage
{

///
/// Holds data for analysis.  Actual analysis is run using the
/// Dataset class.  The Dataset class contains pointers back 
///
class Dataholder
{
public:
  Dataholder();
  ~Dataholder();

  /// Adds an individual to the dataholder
  void add_ind(Individual& ind);
  
  /// Adds an individual to the dataholder
  void add_ind(Individual* ind);

  /// Returns pointer to indiated individual
  inline Individual* operator[] (unsigned int indIndex){return inds[indIndex];}

  /// Returns number of individuals in set
  unsigned int num_inds(){return inds.size();}

  /// Add variable name
  inline void add_covar_name(string var_name){covars.push_back(var_name);covars_map[var_name]=covars.size()-1;}

  /// Retrieve variable name
  inline string get_covar_name(unsigned int index){return covars[index];}

  /// Retrieve index by name
  inline unsigned int get_covar_index(string covar_name){return covars_map[covar_name];}

  /// Number of covariates in set
  inline unsigned int num_covars(){return covars.size();}

  /// Add genotype name
  inline void add_geno_name(string var_name){genos.push_back(var_name);genos_map[var_name]=genos.size()-1;}

  /// Retrieve index for name
  inline unsigned int get_geno_index(string geno_name){return genos_map[geno_name];}

  /// Retrieve variable name
  inline string get_geno_name(unsigned int index){return genos[index];}

  /// Number of genotypes in set
  inline unsigned int num_genotypes(){return genos.size();}
  
  /// Sets maximum locus value in set
  inline void set_max_locus_value(unsigned int max_locus_value){max_locus=max_locus_value;}

  /// Returns maximum  locus value in set
  inline unsigned int get_max_locus_value(){return max_locus;}

  /// Sets indicator when any genotpyes are missing in set
  inline void any_missing_genos(bool missing){any_missing = missing;}

  /// Gets indicator true when mmissing genotypes in data
  inline bool any_missing_genos(){return any_missing;}

  /// Returns genotype for indicated individual at indicated locus
  inline float get_genotype(unsigned int currInd, unsigned int currLoc){return inds[currInd]->value(currLoc);}

  /// Returns number of genotypes in set
  inline unsigned int num_genos(){return inds[0]->num_genotypes();}

  /// Returns number of covariates in set
  inline unsigned int num_covariates(){return inds[0]->num_covariates();}

  /// Returns pointer to ind
  inline Individual* get_ind(unsigned int index){return inds[index];}

  /// Returns individual based on id
  Individual* get_ind_by_id(string id){
    if(inds_map.find(id)==inds_map.end()) throw DataExcept("Unable to find individual with ID=" + id);
    return inds_map[id];
  }
  
  /// Gets value for continous variable missing
  float get_missing_covalue(){return missing_covalue;}
  
  /// Sets missing value for continuous variables
  void set_missing_covalue(float miss){missing_covalue = miss;}
  
  /// Gets missing value for genotypes
  int get_missing_genotype(){return missing_genotype;}
  
  /// Sets missing value for genotypes
  void set_missing_genotype(int miss){missing_genotype = miss;}
  
  /// Adds default snp names
  void add_default_snps();
  
  /// Adds default covariate names
  void add_default_covars();
  
  /// Indicates whether ott-dummy encoding used for genotypes
  void ott_dummy_encoding(bool val){ott_encoded = val;}
  
  /// Returns status on ott_dummy_encoding
  bool ott_dummy_encoding(){return ott_encoded;}
  
  /// Set the cut-off for splitting the training and testing when the 2 files are split
  void set_test_split(int indnum){split_num = indnum;}
  
  /// Returns the split number
  int get_test_split(){return split_num;}
  
private:
  std::vector<Individual*> inds;
  std::vector<std::string> genos;
  std::vector<std::string> covars;
  std::map<std::string, unsigned int> genos_map, covars_map; //key is geno name, value is index into genos array
  std::map<std::string, Individual*> inds_map;
  unsigned int max_locus;
  bool any_missing, ott_encoded;
  int missing_genotype, split_num;
  float missing_covalue;

};

}

#endif /*DATAHOLDER_H_*/
