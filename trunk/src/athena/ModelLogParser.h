//ModelLogParser.h
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

#include <vector>
#include <string>
#include <iostream>
#include "AthenaExcept.h"

struct logModel{
    float fitness;
    std::string model;
    int gen, rank, gram_depth, nn_depth, n_c, n_g;
};

bool compareLogModels(logModel first, logModel second);


///
/// Compiles all Model Log files into a single file per 
/// cross-validation.
///
class ModelLogParser{

    public:
        void compile_files(std::vector<std::string>& filenames, std::string outfilename,
            float not_valid);

    private:
        logModel get_model(std::string& line);
        
        void parse_file(std::string filename, std::vector<std::vector<logModel> >& models);

        void write_output(std::ostream & os, std::vector<std::vector<logModel> >& models,
            float not_valid);

};



