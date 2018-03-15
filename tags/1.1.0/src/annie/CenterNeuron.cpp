#include "Exception.h"
#include "CenterNeuron.h"

#include <sstream>
#include <cmath>
#include <cstdio>
using namespace std;
namespace annie
{

CenterNeuron::CenterNeuron(int label, int dimension) : Neuron(label)
{
	_dimension = dimension;
	_activationFunction = gaussian;
	_dActivationFunction = dgaussian;

	_center = new real[dimension];
	//	_deltaCenter = new real[dimension];
	for (int i=0;i<dimension;i++)
	{
		_center[i] = random();	//random01 ?
		//		_deltaCenter[i]=(real)0.0;
	}

	_classHeirarchy.push_back(_CENTER_NEURON_STRING);
}

CenterNeuron::CenterNeuron(int label, Vector center) : Neuron(label)
{
	_dimension = center.size();
	setCenter(center);
	_activationFunction = gaussian;
	_dActivationFunction = dgaussian;
}

CenterNeuron::CenterNeuron(int label, int dimension, real center[]) : Neuron(label)
{
	_dimension = dimension;
	_activationFunction = gaussian;
	_dActivationFunction = dgaussian;
	setCenter(center);
}

CenterNeuron::~CenterNeuron()
{
	//	delete[] _deltaCenter;
	delete[] _center;
}

void
CenterNeuron::setActivationFunction(ActivationFunction f, ActivationFunction df)
{
	_activationFunction = f;
	_dActivationFunction = df;
}

void
CenterNeuron::setCenter(real center[])
{
	int dim = getDimension();
	for (int i=0;i<dim;i++)
		_center[i]=center[i];
}

void
CenterNeuron::setCenter(Vector center)
{
	if (center.size() != getDimension())
		throw Exception("CenterNeuron::setCenter() - Center provided has different dimension than that of this neuron.");
	Vector::iterator it;
	int i;
	for (i=0, it=center.begin(); it!=center.end(); it++,i++)
		_center[i]=(*it);
}

Vector
CenterNeuron::getCenter() const
{
	Vector center;
	for (uint i=0;i<getDimension();i++)
		center.push_back(_center[i]);
	return center;
}

void
CenterNeuron::_recacheOutput() const
{
	if (_outputCacheValid)
		return;
	if (getInputCount() != getDimension())
		throw Exception("CenterNeuron::_recacheOutput() - Attempt to calculate output even though neuron doesn't have as many inputs as it the center's dimension.");
	_activationCache = 0.0;
	LINKS::const_iterator i;
	int j;
	for (i = _inputLinks.begin(), j = 0; i!=_inputLinks.end(); i++, j++)
	{
		const Link *l = (Link *)(*i);
		real temp = l->getSource()->getOutput() - _center[j];
		_activationCache += (temp*temp);
	}
	_activationCache = (real)sqrt(_activationCache);
	_outputCache = _activationFunction(_activationCache);

	_outputCacheValid = true;
}

void
CenterNeuron::_recacheError() const
{	_errorCacheValid = true;	}

uint
CenterNeuron::getDimension() const
	{	return _dimension;	}

void
CenterNeuron::connect(Neuron *from)
{
	disconnect(from);
	Link *link = new Link(this,from,1.0);
	_inputLinks.push_back(link);
	from->_outputLinks.push_back(link);
	invalidateOutputCache();
}

//void
//CenterNeuron::calculateNewCenter(real learningRate)
//{
//	LINKS::iterator k1,i;
//	int k;
//	real act,dact,factor;
//	Link *inl,*outl;
//	if (getInputCount() != getDimension())
//	{
//		string error(getClassName());
//		error = error + "::calculateNewCenter() - Setup not valid, number of inputs not equal to dimension of center.";
//		throw Exception(error);
//	}
//
//	for (k1=_inputLinks.begin(), k=0; k1!=_inputLinks.end(); k1++,k++)
//	{
//		inl=(Link*)(*k1);
//		_deltaCenter[k] = (real)0.0;
//		act = getActivation();
//		dact = (real)(exp(-1*act*act)*2);//_dActivationFunction(act); //
//		factor = (inl->getSource()->getOutput() - _center[k])*dact;///act;
//		for (i=_outputLinks.begin(); i!=_outputLinks.end(); i++)
//		{
//			outl=(Link*)(*i);
//			_deltaCenter[k] += outl->getDestination()->getError() * outl->getWeight();
//		}
//		_deltaCenter[k] = _deltaCenter[k] * factor * learningRate;
//	}
//}
//
//void
//CenterNeuron::update()
//{
//	int dim = getDimension();
//	int i;
//	for (i=0; i<dim; i++)
//	{
//		_center[i]+=_deltaCenter[i];
//		_deltaCenter[i]=(real)0.0;
//	}
//	invalidateOutputCache();
//	invalidateErrorCache();
//}

const char *
CenterNeuron::getClassName() const
{
	return _CENTER_NEURON_STRING;
}

CenterNeuron::operator string() const
{
	stringstream ans;
	ans << Neuron::operator string();
	ans << "- Center = (";

	for (uint i=0;i<getDimension();i++)
		ans <<_center[i];
	ans << ")\n";
	return ans.str();
}

}; //namespace annie

