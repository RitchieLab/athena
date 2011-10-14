#ifndef OTTDUMMYCONVERT_H_
#define OTTDUMMYCONVERT_H_

#include "Dataholder.h"

namespace data_manage
{


///
/// Converts genotypes in individuals to the ott dummy encoding.
/// Each single genotype is replaced by a pair of genotypes
///
class OttDummyConvert
{
public:
  OttDummyConvert();
  ~OttDummyConvert();

  /// converts all genotypes to the ott-dummy encoding
  void convert_genotypes(Dataholder* holder);

private:
  std::vector<std::vector<char> > convertor;

};

}

#endif /*OTTDUMMYCONVERT_H_*/
