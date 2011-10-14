#include "Exception.h"
#include "Neuron.h"
#include "Link.h"

#include <cmath>
#include <cstdio>
#include <sstream>
using namespace std;

namespace annie
{

real identity(real x)
{	return x;	}

real didentity(real x)
{	return 1.0;	}

real sigmoid(real x)
{
// smd changed to always use the same approximation as in ATHENA
//#ifdef SIGMOID_APPROX	
	if(x < -SIGMOID_APROX_THRESHOLD) return 0;
	if(x > SIGMOID_APROX_THRESHOLD) return 1;
//#endif
	return (real)(1/(1+exp(-x)));
}

real dsigmoid(real x)
{
#ifdef SIGMOID_APPROX	
	if(x < -SIGMOID_APROX_THRESHOLD || x > SIGMOID_APROX_THRESHOLD) return 0;
#endif

	real f=sigmoid(x);
	return f*(1-f);
}

real gaussian(real x)
{
	return (real)(exp(-1*x*x));
}

real dgaussian(real x)
{
	return -2*gaussian(x)*x;
}

real signum(real x)
{
	if (x<0.0)
		return (real)-1.0;
	return (real)1.0;
}

real tansig(real x)
{ return (real)((2/(1+exp(-2*x)))-1); }

real dtansig(real x)
{
	real f=tansig(x);
	return 1-(f*f);
}

Neuron::Neuron(int label)
{
	_label = label;
	_inputLinks.clear();
	_outputLinks.clear();
	_errorCache = _activationCache = _outputCache = 0.0;
	_outputCacheValid = _errorCacheValid = true;
	_activationFunction = identity;
	_classHeirarchy.push_back("Neuron");
}

Neuron::Neuron(Neuron &neuron)
{
	//string error(getClassName());
	string error("Neuron");
	error = error + "::Neuron() - Copy constructor not implemented";
	throw Exception(error);
}

void
Neuron::invalidateOutputCache()
{
	if (_outputCacheValid)
	{
		_outputCacheValid = false;
		LINKS::iterator it;
		for (it=_outputLinks.begin();it!=_outputLinks.end();it++)
		{
			Link *l = (Link *)(*it);
			l->getDestination()->invalidateOutputCache();
		}
	}
}

void
Neuron::invalidateErrorCache()
{
	if (_errorCacheValid)
	{
		_errorCacheValid = false;
		LINKS::iterator it;
		for (it =  _inputLinks.begin(); it != _inputLinks.end(); it++)
		{
			Link *l = (Link *)(*it);
			l->getSource()->invalidateErrorCache();
		}
	}
}

real
Neuron::getActivation() const
{
	_recacheOutput();
	return _activationCache;
}

real
Neuron::getOutput() const
{
	_recacheOutput();
	return _outputCache;
}

real
Neuron::getError() const
{
	_recacheError();
	return _errorCache;
}

uint
Neuron::getInputCount() const
	{	return _inputLinks.size(); }

uint
Neuron::getOutputCount() const
	{	return _outputLinks.size(); }

int
Neuron::getLabel() const
	{	return _label;	}

Neuron::~Neuron()
{
	LINKS::iterator it;
	while (!_inputLinks.empty())
	{
		it = _inputLinks.begin();
		Link *l = (Link *)(*it);
		delete l;
	}
	while (!_outputLinks.empty())
	{
		it = _outputLinks.begin();
		Link *l = (Link *)(*it);
		delete l;
	}
}

Neuron::operator string() const
{
	stringstream ans;
	ans << getClassName() << "(" << getLabel() << ")\n";

	ans << "- Inputs  : " << getInputCount() << "\n";
	if (getInputCount() > 0)
	{
		ans << "  (source,weight)      = ";
		LINKS::const_iterator it;
		for (it=_inputLinks.begin(); it!=_inputLinks.end(); it++)
		{
			Link *l = (Link *)(*it);
			ans << "(" << l->getSource()->getLabel() << "," << l->getWeight() << ") ";
		}
		ans << "\n";
	}

	ans << "- Outputs : " << getOutputCount() << endl;
	if (getOutputCount() > 0)
	{
		ans << "  (destination,weight) = ";
		LINKS::const_iterator it;
		for (it=_outputLinks.begin(); it!=_outputLinks.end(); it++)
		{
			Link *l = (Link *)(*it);
			ans << "(" << l->getDestination()->getLabel() << "," << l->getWeight() << ") ";
		}
		ans << "\n";
	}
	return ans.str();
}

void
Neuron::disconnect(Neuron *from)
{
	//OPT: this is very slow if connect/disconnect is used frequently
	Link *link = new Link(this,from,0.0);
	LINKS::iterator it;
	for (it=_inputLinks.begin();it!=_inputLinks.end();it++)
	{
		Link *l = (Link *)(*it);
		if (l->isEqualTo(link))
		{
			delete l;
			break;
		}
	}
	invalidateOutputCache();
	delete link;
}

int
Neuron::getInputs(vector<int> &labels, Vector &weights)
{
	labels.clear();
	weights.clear();
	LINKS::iterator it;
	for (it = _inputLinks.begin(); it!=_inputLinks.end(); it++)
	{
		Link *l = (Link *)(*it);
		labels.push_back(l->getSource()->getLabel());
		weights.push_back(l->getWeight());
	}
	return getInputCount();
}

void Neuron::getWeights(Vector &out)	const {
    out.resize(_inputLinks.size());
    LINKS::const_iterator it; uint i=0;
    for (it = _inputLinks.begin(); it != _inputLinks.end(); it++, i++) {
    const Link *l = (const Link *)(*it);
	    out[i] = l->getWeight();
    }
}

real
Neuron::getWeight(Neuron *from) const
{
	LINKS::const_iterator it;
	Link *l;
	for (it = _inputLinks.begin(); it!=_inputLinks.end(); it++)
	{
		l = (Link*)(*it);
		if (l->getSource() == from)
			return l->getWeight();
	}
	return 0.0;
}

ostream& operator << (ostream& s, const Neuron &d)
{
	s << (string) d;
	return s;
}

bool
Neuron::instanceOf(const char *className) const
{
	vector<char *>::const_iterator it;
	for (it = _classHeirarchy.begin(); it!= _classHeirarchy.end(); it++)
	{
		char *name = (char*)(*it);
		if (!strcmp(name,className))
			return true;
	}
	return false;
}

} //namespace annie
