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
/*
 * File:   Config.h
 * Author: dudeksm
 *
 * Created on November 5, 2008, 4:59 PM
 */

#ifndef _CONFIG_H
#define	_CONFIG_H

#include <map>
#include <string>
#include <vector>
#include "AthenaExcept.h"
#include "Structs.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_CXX_MPI
#include <mpi.h>
#endif

/// holds parameters and name of the Algorithm
struct AlgorithmParams{
		std::map<std::string, std::string> params;
		std::string name;
};

///
/// Contains configuration parameters for running HEMANN
///
class Config{

public:
		/// Constructor
		Config();

		std::string getDataSetType(){return dataType;}
		void setDataSetType(std::string dtype){dataType = dtype;}

		std::string getOutputName(){return outName;}
		void setOutputName(std::string oname){outName = oname;}

		std::string getMapName(){return mapName;}
		void setMapName(std::string mname){mapName = mname;}

		int getMissingValue(){return missValue;}
		void setMissingValue(int mval, std::string id=""){
			if(mval >= 0 && mval <= 2)
				throw AthenaExcept("The missing value " + id +" cannot be between 0 and 2 inclusive.");
			missValue = mval;
		}

		float getStatusMissingValue(){return statMissValue;}
		void setStatusMissingValue(float val){statMissValue = val;}

		int getRandSeed(){return randSeed;}
		void setRandSeed(int rseed){randSeed = rseed;}

		int getNumCV(){return nCV;}
		void setNumCV(int numCV, std::string id=""){
			if(numCV < 1)
				throw AthenaExcept("Number of cv " + id + " must be greater than zero");
			nCV = numCV;
		}

		std::string getDataSetName(){return dataFile;}
		void setDataSetName(std::string dname){dataFile = dname;}

		void addAlgorithmParam(AlgorithmParams alg){algParams.push_back(alg);}
		std::vector<AlgorithmParams> getAlgorithmParams(){return algParams;}

		std::string getContinFileName(){return continFile;}
		void setContinFileName(std::string cname){continFile = cname;}

		bool getIDinData(){return idIncluded;}
		void setIDinData(bool id){idIncluded = id;}

		float getContinMiss(){return continMiss;}
		void setContinMiss(float val){continMiss = val;}

		bool getOttEncoded(){return ottDummyEncoded;}
		void setOttEncoded(bool ott){ottDummyEncoded = ott;}

		int getNumExchanges(){return numExchanges;}
		void setNumExchanges(int numEx, std::string id=""){
			if(numEx < 0)
				throw AthenaExcept("Number of exchanges " + id + " must be greater than or equal to zero");
			numExchanges = numEx;
		}

		void setCVOutput(bool val){cvOut = val;}
		bool getCVOutput(){return cvOut;}

		void setStatusAdjust(std::string statChange){statusChange=statChange;}
		std::string getStatusAdjust(){return statusChange;}

		void setContinAdjust(std::string cChange){continChange=cChange;}
		std::string getContinAdjust(){return continChange;}

		void setIndOutput(bool io){indsOutput = io;}
		bool getIndOutput(){return indsOutput;}

		void setOutputAllNodesBest(bool val){allNodesOut = val;}
		bool outputAllNodesBest(){return allNodesOut;}

		std::string getTrainFile(){return trainFile;}
		void setTrainFile(std::string filename){trainFile = filename;}

		std::string getTestFile(){return testFile;}
		void setTestFile(std::string filename){testFile = filename;}

		void setContinTestFile(std::string filename){continTest = filename;}
		std::string getContinTestFile(){return continTest;}

		void setContinTrainFile(std::string filename){continTrain = filename;}
		std::string getContinTrainFile(){return continTrain;}

		void setBioFilterFile(std::string filename){bioFilterFile = filename;}
		std::string getBioFilterFile(){return bioFilterFile;}

		void setBioGeneFile(std::string filename){bioGeneFile = filename;}
		std::string getBioGeneFile(){return bioGeneFile;}

		void setBioArchiveFile(std::string filename){bioArchiveFile = filename;}
		std::string getBioArchiveFile(){return bioArchiveFile;}

		void setBioFileType(std::string filetype){biofilterFileType = filetype;}
		std::string getBioFileType(){return biofilterFileType;}

		void setContinMapName(std::string filename){continMap=filename;}
		std::string getContinMapName(){return continMap;}

		void setSelectBestModel(bool tf){selectBest=tf;}
		bool selectBestModel(){return selectBest;}


		void setSplitFile(std::string filename){splitFile = filename;}
		std::string getSplitFile(){return splitFile;}

		void setStartCV(int cv){startCV=cv;}
		int getStartCV(){return startCV;}

		inline std::string getImgWriter(){return imgWriter;}
		inline void setImgWriter(std::string executable){imgWriter=executable;}

		inline std::string getValidationSumFile(){return validationSumFile;}
		inline void setValidationSumFile(std::string filename){validationSumFile=filename;}

		/// throws an exception if parameters are in error
		void checkConfig();

		enum SummaryType{
			True,
			False,
			Best,
			All,
			Suppress
		};

		inline void setSummaryOnly(std::string val){
			std::map<std::string, SummaryType>::iterator iter = summaryMap.find(val);

			if(iter != summaryMap.end()){
				summaryOnly = iter->second;
			}
			else{
				throw AthenaExcept(val + " is not a valid parameter for summary type");
			}
		}

		SummaryType getSummaryOnly(){return summaryOnly;}

		inline void setEncodeType(std::string encodeType){
			encodeName = encodeType;
		}

		inline string getEncodeType(){return encodeName;}


		inline void setLogType(std::string logType){
			std::map<std::string, LogType>::iterator iter = logTypeMap.find(logType);
			if(iter != logTypeMap.end())
				logTypeSelected = iter->second;
			else
				throw AthenaExcept(logType + " is not a valid parameter for log type selection");
		}

		inline LogType getLogType(){return logTypeSelected;}

private:

		void initialize();
		std::map<std::string, SummaryType> summaryMap;
		LogType logTypeSelected;
		std::map<std::string, LogType> logTypeMap;

		std::string dataType, outName, mapName, dataFile, continFile, statusChange,
			trainFile, testFile, continTest, continTrain, bioFilterFile, biofilterFileType,
			bioArchiveFile, bioGeneFile, continMap, splitFile, continChange, encodeName, imgWriter,
			validationSumFile;
		int missValue, nCV, randSeed, numExchanges, startCV;
		float continMiss, statMissValue;
		bool idIncluded, ottDummyEncoded, cvOut, indsOutput, allNodesOut, selectBest;
		std::vector<AlgorithmParams> algParams;
		SummaryType summaryOnly;

};


#endif	/* _CONFIG_H */

