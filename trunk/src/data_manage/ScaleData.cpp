#include "ScaleData.h"
#include <cmath>

namespace data_manage
{

ScaleData::ScaleData()
{
}

ScaleData::~ScaleData()
{
}

///
/// Alters status by putting status through a sigmoid function
/// @param holder Dataholder
///
void ScaleData::sigmoid_status(Dataholder* holder){
  
  Individual * ind; 
  float status;
  
  for(unsigned int curr_ind = 0; curr_ind < holder->num_inds(); curr_ind++){
    ind = (*holder)[curr_ind];
    status = ind->get_status();
    if(status < -709) ind->set_status(-1.0);
    else if(status > 709)ind->set_status(1.0);
    else ind->set_status(1.0/(1.0+exp(-status)));
  }
        
}


}
