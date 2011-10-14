#include "Exception.h"
#include "SimpleNeuron.h"

#include <cstdio>
using namespace std;

namespace annie
{

SimpleNeuron::SimpleNeuron(int label, bool hasBias) : AbstractNeuron(label, hasBias)
{
	_activationFunction = sigmoid;
	_dActivationFunction = dsigmoid;
	_classHeirarchy.push_back(_SIMPLE_NEURON_STRING);
}

void
SimpleNeuron::setActivationFunction(ActivationFunction f, ActivationFunction df)
{
	_activationFunction = f;
	_dActivationFunction = df;
	invalidateOutputCache();
	invalidateErrorCache();
}

void
SimpleNeuron::setDesiredOutput(real desired)
{
	if (getOutputCount()!=0)
	{
		string error(getClassName());
		error = error + "::setDesiredOutput() - Called for a non-output neuron";
		throw Exception(error);
	}

	invalidateErrorCache();
	// \delta_j = (d_j - y_j) * dy_j
	_errorCache =  (desired - getOutput()) * _dActivationFunction(getActivation());
	VM("delta_" << getLabel() << " = " << "(" << desired << "-" <<  getOutput() << ") * " << "f'(" << getActivation() << ")" << " = " << _errorCache << endl)
	_recacheError();	//this only sets _errorCacheValid = true, because this is an output neuron
}

void
SimpleNeuron::_recacheError() const
{
	if (getOutputCount()==0) //ie an output neuron
		_errorCacheValid = true;
	if (_errorCacheValid)
		return;
	//Note that this work will be done only if this is a NON-OUTPUT neuron
	//ie, it's output is taken as input by other neurons

	//\delta_j = (\sum_k ( \delta_k w_jk)) * dy_j

	LINKS::const_iterator it;
	_errorCache = 0;
	for (it = _outputLinks.begin(); it!=_outputLinks.end(); it++)
	{
		const Link *l = (Link *)(*it);
		_errorCache += l->getWeight() * l->getDestination()->getError();
	}
	
	_errorCache *= _dActivationFunction(getActivation());
	_errorCacheValid = true;

#ifdef VM_ENABLED
	VM("delta_" << getLabel() << " = ( ")
	//real ec = 0;
	for (it = _outputLinks.begin(); it!=_outputLinks.end(); it++)
	{
		const Link *l = (Link *)(*it);
		if(it != _outputLinks.begin()) {VM(" + ");}
		VM(	l->getWeight() << " * " << l->getDestination()->getError());
		//ec += l->getWeight() * l->getDestination()->getError();
	}
	VM(" ) * f'( " << getActivation() << ") = " << getError() << endl);
#endif
}

void
SimpleNeuron::_recacheOutput() const
{
	if (_outputCacheValid)
		return;
	LINKS::const_iterator it;

	_activationCache = 0;

	for (it=_inputLinks.begin();it!=_inputLinks.end();it++)
	{
		const Link *l = (Link *)(*it);
		_activationCache += l->getSource()->getOutput() * l->getWeight();
	}
	
	if (_hasBias) _activationCache += _bias;
	_outputCache = _activationFunction(_activationCache);
	_outputCacheValid = true;

#ifdef VM_ENABLED
	//must make special run ~~ so that the outputs from individual neurons don't mix..
	VM("y_" << getLabel() << " = f(")
	for (it=_inputLinks.begin();it!=_inputLinks.end();it++)
	{
		const Link *l = (Link *)(*it);
		VM(l->getSource()->getOutput() << "*" <<  l->getWeight())
		VM(" + ")
	}
	if (_hasBias) VM(_bias);

	VM(") = f(" << _activationCache << ") = " << _outputCache << endl)
#endif	//VM
}

const char *
SimpleNeuron::getClassName() const
{
	return _SIMPLE_NEURON_STRING;
}

} //namespace annie


