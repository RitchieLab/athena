// BestModelSelector.h

#ifndef _BIOMODELSELECTOR_H
#define	_BIOMODELSELECTOR_H


#include<vector>
#include<map>
#include<set>
#include<iostream>
#include <Spearman.h>
#include "Solution.h"
#include<Utility.h>


class BestModelSelector{
	public:
	
	BestModelSelector();
	
	/// find best variables from the solutions passed
	void selectBestVariables(std::vector<Solution*> & solutions, Dataholder * data);

	friend ostream& operator<< (ostream& os, BestModelSelector & selector);
	
	inline void setCorrThreshold(float val){correlationThreshold = val;}
	
	inline float getCorrThreshold(){return correlationThreshold;}
	
	inline void setCVThreshold(int cv){crossValThreshold=cv;}
	
	inline std::vector<int> getIncludedGenos(){return includeGenos;}
	
	inline std::vector<int> getIncludedContins(){return includeContins;}
	
	struct CorrelationScore{
			int var1, var2;
			std::string name1, name2;
			float score;
		};
	
  inline std::vector<CorrelationScore>& getGenoCorrelations(){return genoCorrelations;}
  
  inline std::vector<CorrelationScore>& getContinCorrelations(){return continCorrelations;}
	
	private:
	
		std::set<int> countVariables(std::vector<int>& vars, 
			std::map<int, int>& varTotals);
			
		void addSignalCV(std::set<int>& varMult,
			std::map<int, int>& varTotals, map<int, vector<int> >& correlationSignals);
		
		std::map<int, vector<int> > calculateCorrelation(std::vector<int>& vars, 
			Dataholder* data, bool continUsed=false);

		void fillGenoPairs(std::vector<stat::PairedValues>& pairs, Dataholder& data, int var1,
			int var2);
			
		void fillContinPairs(std::vector<stat::PairedValues>& pairs, Dataholder& data, int var1,
			int var2);	
	

	
		float correlationThreshold;
		int crossValThreshold;
		std::vector<int> includeContins, includeGenos;
		std::vector<CorrelationScore> genoCorrelations, continCorrelations;
};


#endif

