/**
 * Functions relating creation of time-series examples
 * @author OP
 */

#include "examples.h"
#include "Exception.h"
#include "TrainingSet.h"
#include "random.h"
#include "Control.h"
#include <cmath>

namespace annie	{

Vector toPolar(const Vector &in)	{
	Vector out(in.size());
	switch(in.size())	{
	case 2:	{
			real r =  in[0];
			real a =  in[1];
			out[0] = 0.5 + 0.5 * r * cos(a * 2. * M_PI);
			out[1] = 0.5 + 0.5 * r * sin(a * 2. * M_PI);
		} break;
	default:
		throw Exception("TODO: heck.. only implemented for 2 dimensions :(");
	}
	return out;
}

Vector normalize(const Vector &in)	{
	Vector r = in;
	r.normalize();
	return r;
}

extern Vector vectorIdentity(const Vector &in) { return in; }

TrainingSet *transformInputs(TrainingSet &source, XformFunction f, uint res_ins)	throw (){
	TrainingSet *ret = new TrainingSet(res_ins, source.getOutputSize());
	Vector in(res_ins), out(source.getOutputSize());
	source.initialize();
	while(!source.epochOver())	{
		source.getNextPair(in, out);
		ret->addIOpair(f(in), out);
	}
	return ret;
}


TrainingSet *uniformCube(uint samples, uint dim)	{
	return randomSamples(samples, dim, random01);
}

TrainingSet *uniformSphere(uint samples, uint dim)	{
	/*
	* NOTE: as the interval is [0, 1] and not [0, 1), this is not exactly uniform ..*/
	TrainingSet *tmp = uniformCube(samples, dim);
	TrainingSet *ret;
	try	{
		ret=transformInputs(*tmp, toPolar);
	} catch(...)	{ delete(tmp); throw; }
	delete(tmp);
	return ret;
}

TrainingSet *randomSamples(uint samples, uint dim, SamplingFunction samplingFunction)	{
	SamplingFunction *sfs = new SamplingFunction[dim];
	for(uint i=0; i<dim; i++) sfs[i] = samplingFunction;
	TrainingSet *ret;
try	{ ret = randomSamples(samples, dim, sfs);	}
	catch(...)	{ delete[] sfs; throw; }
	delete[] sfs;
	return ret;
}

TrainingSet *randomSamples(uint samples, uint dim, SamplingFunction *samplingFunctions)	{
	TrainingSet *res = new TrainingSet(dim, 0);
	Vector in(dim), out(0);
	for(uint i=0; i<samples; i++)	{
		for(uint j=0; j<dim; j++) in[j] = samplingFunctions[j]();
		res->addIOpair(in, out);
	}
	return res;
}

TrainingSet *transformInputs(TrainingSet &source, XformFunction f) throw()	{
	return transformInputs(source, f, source.getInputSize());
}

void Selector::xform(const Vector &in1, const Vector &out1, Vector &in2, Vector &out2) const	{
	for(uint i=0; i<_ins.size(); i++)	{
		int s = _ins[i];
		in2[i] = s >= 0 ? in1[s] : out1[-s - 1];
	}
	
	for(uint i=0; i<_outs.size(); i++)	{
		int s = _outs[i];
		out2[i] = s >= 0 ? in1[s] : out1[-s - 1];
	}
}

}	//annie
