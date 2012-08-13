//ModelLogParser.h

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
        void compile_files(std::vector<std::string>& filenames, std::string outfilename);

    private:
        logModel get_model(std::string& line);
        
        void parse_file(std::string filename, std::vector<std::vector<logModel> >& models);

        void write_output(std::ostream & os, std::vector<std::vector<logModel> >& models);

};



