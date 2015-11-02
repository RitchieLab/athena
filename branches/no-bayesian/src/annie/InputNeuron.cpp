#include "InputNeuron.h"
#include <sstream>

using namespace std;
namespace annie
{

InputNeuron::InputNeuron(int label) : Neuron(label)
{
	_classHeirarchy.push_back(_INPUT_NEURON_STRING);
}

void
InputNeuron::_recacheOutput() const
{	_outputCacheValid = true;	}

void
InputNeuron::_recacheError() const
{	_errorCacheValid = true;	}

void
InputNeuron::setValue(real value)
{
	_activationCache = _outputCache = value;
	invalidateOutputCache();
	_recacheOutput();
}


InputNeuron::operator string() const
{
	stringstream out;
	out << "Input" << (string)(Neuron&)(*this) << ", value: " << _outputCache;
	return out.str();
}

const char *
InputNeuron::getClassName() const
{
	return _INPUT_NEURON_STRING;
}

}; //namespace annie

