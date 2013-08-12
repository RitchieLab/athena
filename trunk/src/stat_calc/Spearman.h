#ifndef SPEARMAN_H
#define SPEARMAN_H

#include <vector>

namespace stat
{

struct PairedValues{
	float values[2];
	float ranks[2];
	
	void setPair(float x, float y){
		values[0] = x;
		values[1] = y;
	}
};


struct RankValues
{
	RankValues(int i){index = i;}

	inline bool operator()(const PairedValues& lv, const PairedValues& rv)
	{
		return lv.values[index] < rv.values[index];
	}
	
	int index;
};


class Spearman
{

private:

	Spearman();
  ~Spearman();

public:

  /// returns the Spearman Rank coefficient for the values passed
  static float calculate(std::vector<PairedValues>& pairs);
  
  /// adds ranks to the pairs
  static void rank(std::vector<PairedValues>& pairs, int index);

};

}

#endif /*SPEARMAN_H*/