#ifndef _ABSTRACTNEURON_H
#define _ABSTRACTNEURON_H

#include "Neuron.h"
#include "defines.h"

namespace annie
{

/**
 * Implementation helper for common types of Neurons. Has everything except activation function,
 * which is defined in offsprings.
 *
 * Note: rather than adding more virtual functions for activation function and derivation, we make virtual
 * recaching functions. This should be more efficient in the performance-critical parts.
 */
class AbstractNeuron : public Neuron
{
protected:
	/// Is this neuron allowed to have a bias?
	bool _hasBias;

	/// If allowed to have a bias then the bias, otherwise 0.0
	real _bias;

	/// The change in bias, calculated using the gradient descent rule
	real _deltaBias;

	/// @see Neuron::_recacheOutput
	virtual void _recacheOutput() const =0;

	/// @see Neuron::_recacheError
	virtual void _recacheError() const =0;

	static const real INIT_WEIGHT_MAX;	//< initial weights will be uniformly sampled from [-INIT_WEIGHT_MAX,  INIT_WEIGHT_MAX]
public:

	static real getRandomWeight();
	/** Creates a simple neuron with the given label.
	  * @param label The label to be given to the neuron
	  * @param hasBias true if the neuron is allowed to have a bias, false otherwise. Default is true
	  * @see removeBias
	  */
	AbstractNeuron(int label, bool hasBias = true);

	/** Sets the bias of the neuron.
	  * @param bias The bias to be given to the neuron
	  * @throws Exception if the neuron is not allowed to have a bias
	  */
	virtual void setBias(real bias);

	/// Is the neuron allowed to have a bias?
	virtual bool hasBias() const;

	/// The bias of the neuron, 0.0 if it's not allowed to have a bias
	virtual real getBias() const;

	/// Sets bias to 0.0 and prevents the neuron from having a bias
	virtual void removeBias();

	/** Sets the desired output of the neuron.
	  * Should be called only for output neurons, i.e., those whose output is not
	  * connected to anyone else. Setting the desired output at these neurons
	  * will form the basis of error backpropagation
	  * @param desired The desired output of this neuron
	  * @throws Exception if the neuron is not an output neuron
	  */
	virtual void setDesiredOutput(real desired)=0;

	/** Connects the given neuron to this one, i.e., the output of the supplied neuron
	  * will be given as input to this neuron. A random weight is provided to the link
	  * @param from The neuron whose output is to be taken is as input
	  */
	virtual void connect(Neuron *from);
	void connect(Neuron &from) { connect(&from); }
	void connect(Neuron &from, real w) { connect(&from, w); }

	/// TODO: rename the above...
	void connectFrom(Neuron *from)	{ connect(from); }

	/** Connects the given neuron to this one, i.e., the output of the supplied neuron
	  * will be given as input to this neuron. The weight of the link will be the one
	  * supplied.
	  * @param from The neuron whose output is to be taken as input
	  * @param weight The weight of the connection
	  */
	virtual void connect(Neuron *from, real weight);

	/** Calculates the adjustment to incoming weights based on the gradient
	  * descent rule (backpropagation).
	  * Doesn't actually update the weights, just sets the value.
	  * @see update
	  */
	virtual void calculateNewWeights(real learningRate, real momentum);

	/** Updates the weights of incoming connections according to the values
	  * calculated using the gradient descent rule.
	  * @see calculateNewWeights
	  */
	virtual void update();

	/// Returns "AbstractNeuron"
	virtual const char *getClassName() const;

	/** Returns the weight of the link between the provided neuron and this neuron.
	  * @param from The neuron whose output is connected to this neuron's input
	  * @return The weight of the connection, 0.0 if no connection exists
	  */
	virtual real getWeight(Neuron *from) const;

	void randomizeWeights();

	/// @see Neuron::string()
	virtual operator std::string() const;
};

}; //namespace annie
#endif // define _ABSTRACTNEURON_H

