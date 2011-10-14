#ifndef _STEPHENDUMMYCONVERT_H_
#define _STEPHENDUMMYCONVERT_H_

#include "Dataholder.h"

namespace data_manage
{


///
/// Converts genotypes in individuals to the ott dummy encoding.
/// Each single genotype is replaced by a pair of genotypes
///
class StephenDummyConvert
{
public:
  StephenDummyConvert();
  ~StephenDummyConvert();

  /// converts all genotypes to the ott-dummy encoding
  void convert_genotypes(Dataholder* holder);

private:
//   std::vector<std::vector<char> > convertor;

};

}

#endif /*_STEPHENDUMMYCONVERT_H_*/
