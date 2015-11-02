#include "Layer.h"
#include "Exception.h"
#include "AbstractNeuron.h"

using namespace std;

namespace annie
{

const int Layer::MAX_LAYER_SIZE = MAX_NEURONS_IN_LAYER;

Layer::Layer(int label)
{
	_size = 0;
	_label = label;
}

Layer::~Layer()
{
	vector<Neuron *>::iterator it;
	while (!_neurons.empty())
	{
		it = _neurons.begin();
		delete *it;
		it = _neurons.erase(it);
	}
}

Vector
Layer::getActivation()
{
	Vector answer;
	vector<Neuron *>::iterator it;
	for (it = _neurons.begin(); it!=_neurons.end(); it++)
	{
		Neuron *n = (Neuron *)(*it);
		answer.push_back(n->getActivation());
	}
	return answer;
}

Vector
Layer::getOutput()
{
	Vector answer;
	vector<Neuron *>::iterator it;
	for (it = _neurons.begin(); it!=_neurons.end(); it++)
	{
		Neuron *n = (Neuron *)(*it);
		answer.push_back(n->getOutput());
	}
	return answer;
}

Neuron &
Layer::getNeuron(uint i)
{
	if (i<0 || i>=getSize())
	{
		string error(getClassName());
		error = error + "::getNeuron() - Invalid index specified";
		throw Exception(error);
	}
	return *(_neurons[i]);
}


const Neuron &
Layer::getNeuron(uint i) const	{
	if (i<0 || i>=getSize())
	{
		string error("zz"); ///XXX
		error = error + "::getNeuron() - Invalid index specified";
		throw Exception(error);
	}
	return *(_neurons[i]);

}//XXX

uint
Layer::getSize() const
	{	return _size;	}

int
Layer::getLabel() const
	{	return _label;	}

void
Layer::addNeuron(Neuron *nrn)
{
	_neurons.push_back(nrn);
	_size++;
}

const char *
Layer::getClassName()
{	return "Layer";	}

}; //namespace annie

