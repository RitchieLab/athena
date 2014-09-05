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
#include "Config.h"

using namespace std;

/// /// Constructor ///
Config::Config(){
		initialize();
}


/// /// Initializes variables ///
void Config::initialize(){
		dataType = "TEXT";
		outName = "athena";
		mapName = "";
		statusChange = "";
		randSeed = 1;
		dataFile = "";
		idIncluded = false;
		continFile = "";
		continMiss = -9999;
		ottDummyEncoded = false;
		numExchanges = 2;
		cvOut = false;
		statusChange = "NONE";
		continChange = "NONE";
		indsOutput = false;
		allNodesOut = false;
		testFile = "";
		trainFile = "";
		continTest = "";
		continTrain = "";
		continMap = "";
		missValue = -1;
		bioFilterFile = "";
		summaryOnly = False;
		statMissValue = -1.0;
		biofilterFileType = "TEXT";
		bioArchiveFile = bioGeneFile = "";
		encodeName = "NONE";
		logTypeMap["NONE"] = LogNone;
		logTypeMap["SUMMARY"] = LogSummary;
		// changed to match Summary as part of log optimization
		logTypeMap["DETAILED"] = LogDetailed;
		logTypeMap["VARIABLES"] = LogVariables;
		logTypeMap["OVERVIEW"] = LogOverview;
		logTypeSelected = LogNone;
		summaryMap["TRUE"] = True;
		summaryMap["FALSE"] = False;
		summaryMap["BEST"] = Best;
		summaryMap["ALL"] = All;
		imgWriter="";
		startCV = 1;
		selectBest=false;
}

/// /// Checks configuration parameters for errors /// @throws AthenaExcept /// void 
void Config::checkConfig(){
	if(numExchanges < 0)
		throw AthenaExcept("NUMEXCHANGES must be greater than zero");
}
