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
 * File:   ConfigFileReader.h
 * Author: dudeksm
 *
 * Created on November 5, 2008, 4:58 PM
 */

#ifndef _CONFIGFILEREADER_H
#define	_CONFIGFILEREADER_H

///
/// Reads configuration file and sets config 
/// objects
///

#include "Config.h"
#include "AthenaExcept.h"
#include <Stringmanip.h>
#include <fstream>

using namespace data_manage;

class ConfigFileReader{
public:
    ConfigFileReader(){initialize_keywords();}
    
    ConfigFileReader(string configfile);
    
    Config read_config(string configfile);
    
private:
    
    /// Reads parameters for current algorithm and stores as map in Config
    void read_algo_params(std::map<std::string, std::string> algo_params, std::ifstream & configFile);
    
    void initialize_keywords();
    
    bool skip_line(std::string line);
    
    void read_alg_params(AlgorithmParams& alg_param, std::ifstream& configstream);
    
    bool param_true(std::string param);
    
    enum configKeyWords{
        keyNoMatch,
        keyDataset,
        keyOut,
        keyMapFile,
        keyAlgorithm,
        keyEnd,
        keyMissingValue,
        keyRandSeed,
        keyDatasetType,
        keyCV,
        keyContinFile,
        keyIDIncluded,
        keyContinMiss,
        keyDummyEncode,
        keyNumExchanges,
        keyCVOutput,
        keyStatusChange,
        keyBestModelIndOutput,
        keyOutputAllNodeBest,
        keyTestFile,
        keyTrainFile,
        keyContinTestFile,
        keyContinTrainFile,
        keyBioFilterFile,
        keySummaryOnly,
        keyStatusMissingValue,
        keyBioFileType,
        keyBioArchiveFile,
        keyBioGeneFile,
        keyLogType,
        keyContinMapFile,
        keySplitFile,
        keyCVStart
    };
    
    std::map<std::string, configKeyWords> keywordMap;
    
    std::string output_name, dataset_type;
    
};

#endif	/* _CONFIGFILEREADER_H */

