#include "Exception.h"

using namespace std;

namespace annie
{

Exception::Exception(const char *info)
{	details.assign(info);	}

Exception::Exception(string info)
{	details=info;	}

}; //namespace annie

