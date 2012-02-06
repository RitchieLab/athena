#include "Exception.h"
#include "AbstractNeuron.h"
#include "random.h"
#include "auxf.h"

#include <cstdio>
#include <sstream>
#include <iostream>
using namespace std;

namespace annie
{

const real AbstractNeuron::INIT_WEIGHT_MAX = 0.5;
AbstractNeuron::AbstractNeuron(int label, bool hasBias) : Neuron(label)
{
	_hasBias = hasBias;
	_bias = _deltaBias = 0.0;
	_classHeirarchy.push_back(_ABSTRACT_NEURON_STRING);
}

void
AbstractNeuron::connect(Neuron *from, real weight)
{
	Link *link = new Link(this,from,weight);
	disconnect(from);
	_inputLinks.push_back(link);
	from->_outputLinks.push_back(link);
	invalidateOutputCache();
}

void
AbstractNeuron::connect(Neuron *from)
{
	connect(from, getRandomWeight());
}

real AbstractNeuron::getRandomWeight()	{
	return uniformRandom(-INIT_WEIGHT_MAX, INIT_WEIGHT_MAX);
}

void
AbstractNeuron::calculateNewWeights(real learningRate, real momentum)
{
	if (getInputCount()==0)
	{
		string error(getClassName());
		error += "::calculateNewWeights() - Called for a neuron with no inputs";
		throw Exception(error);
	}
	LINKS::iterator it;
	//VM("deltas for neuron " << getLabel() << endl)
	for (it=_inputLinks.begin();it!=_inputLinks.end();it++)
	{
		Link *link = *it;
		real input = link->getSource()->getOutput();
		/*
		 * \delta(w_ij)(t) = learningRate * \delta_j * y_j + momentum * (w_ij(t) - w_ij(t-1))
		 */
		link->setDeltaWeight(getError() * input * learningRate  +  momentum * link->getLastDeltaWeight());
		VM("\tDw_" << getLabel() << "," << link->getSource()->getLabel() << " = " << getError() << " * " <<  input << " * " <<  learningRate << " = " << getError() * input * learningRate << endl)
		if(momentum != 0)
			VM("\t .. + " << momentum << " * " << link->getLastDeltaWeight() << " = " << getError() * input * learningRate  +  momentum * link->getLastDeltaWeight());
	}
	if (_hasBias)	{
		_deltaBias = getError()*learningRate;	//TODO: momentum not applied to bias ..
		VM("Dbias_" << getLabel() << " = " << getError() << " * " << learningRate << " = " << _deltaBias);
	}
	VM(endl);
}

void
AbstractNeuron::update()
{
	LINKS::iterator it;
	for (it=_inputLinks.begin();it!=_inputLinks.end();it++)
	{
		Link *link = (Link*)(*it);
		link->updateWeight();
	}
	if (_hasBias)
	{
		_bias+=_deltaBias;
		VM("\tbias_" << getLabel() << " = " << _bias << "\n");
		_deltaBias=0.0;
	}
	invalidateOutputCache();
	invalidateErrorCache();
}

void
AbstractNeuron::removeBias()
{
	_hasBias = false;
	_bias = 0.0;
	invalidateOutputCache();
	invalidateErrorCache();
}

void
AbstractNeuron::setBias(real bias)
{
	if (!_hasBias)
		throw Exception("AbstractNeuron::setBias() - This neuron isn't supposed to have a bias");
	_bias = bias;
	invalidateOutputCache();
}

real
AbstractNeuron::getBias() const
{	return _bias;	}

bool
AbstractNeuron::hasBias() const
{	return _hasBias;	}

AbstractNeuron::operator string() const {
	stringstream ans;
	ans << Neuron::operator string();
	if (_hasBias) ans << "- Bias: " << getBias() << endl;
	return ans.str();
}

const char *
AbstractNeuron::getClassName() const
{
	return _ABSTRACT_NEURON_STRING;
}

real
AbstractNeuron::getWeight(Neuron *from) const
{
	LINKS::const_iterator i;
	for (i=_inputLinks.begin(); i!=_inputLinks.end(); i++)
	{
		Link *l = (Link*)(*i);
		if (l->getSource() == from)
			return l->getWeight();
	}
	return 0.0;
}

void AbstractNeuron::randomizeWeights()	{
	LINKS::iterator it;
	for (it = _inputLinks.begin(); it!=_inputLinks.end(); it++) {
		Link *l = (Link *)(*it);
		l->setWeight(getRandomWeight());
	}

	if(hasBias()) setBias(getRandomWeight());
}

} //namespace annie


