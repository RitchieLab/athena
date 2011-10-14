#ifndef E_T_H
#define E_T_H

#include <vector>
#include "defines.h"
#include "TrainingSet.h"

namespace annie	{

/**
 * Function that can be sampled in discrete time steps
 */
typedef real (*TimeSeriesCreatingFunction)(unsigned time);

/**
 * Store for samples
 */
typedef std::vector<real> SamplesContainer;


/**
 * Sample the given function
 */
void sampleFunction(TimeSeriesCreatingFunction f, unsigned from, unsigned to, SamplesContainer &out);

/**
 * Make time-series examples from the given sample set.
 *
 * Examples will look like s[j],...,s[j+inputs-1] --> s[j+inputs],...,s[j+inputs+outputs-1]
 */
class TrainingSet;
TrainingSet *makeExamples(unsigned inputs, unsigned outputs, const SamplesContainer &samples, int from=0, int to=-1);


/// generates one sample at each call
/// see... we have random01 :)
typedef real (*SamplingFunction)();


/**
 * Produce set of points randomly drawn from the N-d cube.
 * Each I-O pair will only have the "I" component - dimension of output will be 0.
 * Each dimension is supposed to be a [0, 1] interval
 * This function is useful for Kohonen Maps, for example
 *
 * @param samples - how many samples do you want?
 * @param samplingFunctions	- array [dim] of sample generators for each dimension
 * @return a newly created training set. Caller is responsible for deletion!!
 */
TrainingSet *randomSamples(uint samples, uint dim, SamplingFunction *samplingFunctions);

/// sugar for randomSamples - apply the same sampling function to all dimensions
TrainingSet *randomSamples(uint samples, uint dim, SamplingFunction samplingFunction);

//TODO: somewhat outdated by TS functionalites ..
/**
 *	Transform the input component of training examples by given function
 * */
TrainingSet *transformInputs(TrainingSet &source, XformFunction f, uint res_ins) throw();

///sugar
TrainingSet *transformInputs(TrainingSet &source, XformFunction f) throw();
/*
 * is this version needed somewhere?
 * transform one TrainingSet to another on sample-by-sample basis
 * 
 * @param res_ins, res_outs specify the dimensions of output. If not specified, the same as source is assumed
 */
//TrainingSet *transformExamples(const TrainingSet &source, XformFunction f, Xform.., uint res_ins=in.getInputCount(), uint res_outs=in.getInputCount());

/**
 *	some basic sampling and transforming functions
 */


/**
 * map input to sphere
 * [0] = r
 * [1] = angle ( <- [0, 1], not radians!!!)
 *
 * actually, the output is scaled to [0, 1]  ( .5 + .5 * )
 */

extern Vector toPolar(const Vector &in);
extern Vector normalize(const Vector &in);
extern Vector vectorIdentity(const Vector &in);

/**
 * Heck, float cannot be used as a template parameter...  probably to decline the piggy things like this one
 * we should use xforming class instead of xforming function (like trhat for mixed i-o's)
 */
template<int T_multiplier1000>
static Vector Xscale(const Vector &in)	{
	return in * ((real)T_multiplier1000 / 1000);
}

template<int T_add1000>
static Vector Xadd(const Vector &in)	{
	return in + ((real)T_add1000 / 1000);
}

/**
 * some shortcuts
 * these are good examples of samplers usage
 */
TrainingSet *uniformCube(uint samples, uint dim);

/// each spike is equally covered => center is denser
TrainingSet *uniformSphere(uint  samples, uint dim);

/**
 * Select features from the TS vectors
 */
struct Selector	: TSTransformer {
	
	/**
	 * At any place of selection, there can be number from the interval [- sourceTS.getOuputCount(); sourceTS.getInputCount()].
	 * Examples:
	 * First vector from output: -1  second: -2
	 * First vector from input: 0 second: 1
	 * The results will be drawn from these features in given order. One index may be specified in this order.
	 */	
	typedef std::vector<int> Selection;
	Selector(Selection ins, Selection outs) : TSTransformer(ins.size(), outs.size()), _ins(ins), _outs(outs)	{}	//TODO: check the bounds
	virtual void xform(const Vector &in1, const Vector &in1, Vector &out1, Vector &out2) const;
  protected:
	Selection _ins, _outs;
};


struct Shrinker : Selector	{
	Shrinker(uint inputs, uint outputs=0) : Selector(interval(0, inputs), interval(-1, -outputs - 1)) {}
  protected:
	/**
	 * @param to not included
	 * @param from included
	 * to < 0 <=> from < 0 !
	 */
	static Selection interval(int from, int to)	{
		Selection res;
		if(from < 0) {
			ASSERT(to <= 0 && from >= to);
			for(int i=from; i>to; i--) res.push_back(i);
		}	else	{
			ASSERT(from <= to);
			for(int i=from; i<to; i++) res.push_back(i);
		}
		return res;
	}

};

} //namespace annie
#endif //_H
