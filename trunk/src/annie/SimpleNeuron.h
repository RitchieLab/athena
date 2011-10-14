#ifndef _SIMPLENEURON_H
#define _SIMPLENEURON_H

#include "AbstractNeuron.h"
#include "defines.h"

namespace annie
{

/** A simple perceptron - i.e., it takes as input the weighted sum of the
  * outputs of the neurons connected to it.
  * Probably the first type of neuron you would come across in any 
  * introductory literature of neural networks.
  */
class SimpleNeuron : public AbstractNeuron
{
protected:
	/// The derivative of the activation function, required for gradient descent training
	ActivationFunction _dActivationFunction;

	/// @see Neuron::_recacheOutput
	virtual void _recacheOutput() const;

	/// @see Neuron::_recacheError
	virtual void _recacheError() const;

public:
	/** Creates a simple neuron with the given label.
	  * @param label The label to be given to the neuron
	  * @param hasBias true if the neuron is allowed to have a bias, false otherwise. Default is true
	  * @see removeBias
	  */
	SimpleNeuron(int label, bool hasBias = true);

	/** Sets the desired output of the neuron.
	  * Should be called only for output neurons, i.e., those whose output is not
	  * connected to anyone else. Setting the desired output at these neurons
	  * will form the basis of error backpropagation
	  * @param desired The desired output of this neuron
	  * @throws Exception if the neuron is not an output neuron
	  */
	virtual void setDesiredOutput(real desired);

	/** Sets the activation function and its derivative (required for error backpropagation)
	  * @param f The activation function to be used
	  * @param df The derivative of the activation function
	  */
	virtual void setActivationFunction(ActivationFunction f, ActivationFunction df);

	/// Returns "SimpleNeuron"
	virtual const char *getClassName() const;
};

}; //namespace annie
#endif // define _SIMPLENEURON_H

