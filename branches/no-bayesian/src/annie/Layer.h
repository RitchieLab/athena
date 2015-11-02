#ifndef _LAYER_H
#define _LAYER_H

#include "Neuron.h"
#include "Link.h"

#include <vector>

namespace annie
{
/** Abstraction for a "layer" of neurons, i.e., a group of neurons not
  * connected to each other.
  * @see InputLayer
  */
class Layer
{
protected:
	/// The label of the layer
	int _label;

	/** The number of neurons in the layer.
	  * If you create a sub-class of this class, then the onus
	  * of ensuring that this value is consistent lies on you!
	  */
	uint _size;

	/// The neurons in this layer.
	std::vector<Neuron *> _neurons;

public:
	/** The maximum number of neurons in a layer
	  * Needed for some automatic label assignments of neurons and layers
	  * in Networks
	  */
	static const int MAX_LAYER_SIZE;

	/** Constructs a layer with the given label */
	Layer(int label);

	///deletes all neurons
	virtual ~Layer();

	/// Returns the label of the layer
	virtual int getLabel() const;

	/// The size of the layer (number of neurons in it)
	virtual uint getSize() const;

	/// complies w/ STL ..
	uint size() const { return getSize(); }

	/// Adds the given neuron to the layer
	/// Layer is responsible for deletion
	virtual void addNeuron(Neuron *nrn);

	/** Gives the ith reference in the layer.
	  * @param i The index of the neuron in the layer (0<=i<getSize()).
	  * @throws Exception if the index given is invalid
	  */
	virtual Neuron& getNeuron(uint i);
	virtual const Neuron& getNeuron(uint i) const;

	Neuron& operator[](uint i)	{ return getNeuron(i); }
	const Neuron& operator[](uint i)	const { return getNeuron(i); }

	/// The activation vector formed by the activations of individual neurons in the layer
	virtual Vector getActivation();

	/// The output vector formed by the outputs of individual neurons in the layer
	virtual Vector getOutput();

	/// Returns "Layer"
	virtual const char *getClassName();
};

///nice getters...
///TODO: typechecking || forcing in addNeurons
template <class T>
class TLayer : public Layer	{
  public:
	TLayer(int label) : Layer(label) {}
	T& operator[](uint i)	{ return (T &) getNeuron(i); }
	const T& operator[](uint i)	const { return (const T &) getNeuron(i); }
};

}; //namespace annie
#endif // define _LAYER_H

