#ifndef MULTILAYERNETWORK_H
#define MULTILAYERNETWORK_H


#include "Layer.h"
#include "InputLayer.h"
#include "Network.h"
#include "TrainingSet.h"
#include "Control.h"

namespace annie
{
/** Abstraction of a multi-layer perceptron network.
  * Basically, layers of SimpleNeurons constitute this network.
  * Training is done using the backpropagation technique which uses
  * the gradient descent method.
  *
  * The labels of the layers start from 0 (for the input layer) and then
  * keep moving on. The labels of neurons in the layers is =
  * <layer label>*Layer::MAX_LAYER_SIZE + <neuron index in layer>
  *
  * All neurons in the layers are allowed to have a bias.
  * 
  * @todo The copy constructor
  */
class MultiLayerNetwork : public Network
{
protected:
	/** Number of layers in the network.
	  * If you derive your own network from this class, the onus of
	  * keeping _nLayers consistent lies on you!
	  * Doesn't count the input "layer"
	  */
	uint _nLayers;

	/// The layers
	std::vector<Layer *> _layers;

	/// The input layer
	InputLayer *_inputLayer();

	/// The output layer
	Layer *_outputLayer();

	/// impl detail - throws then layer is not valid
	void _layerValid(uint layer) const;

	void _connectLayer(Layer &srcLayer, Layer &destLayer);
public:
	static const real DEFAULT_MOMENTUM, DEFAULT_LEARNINGRATE;
	static const Creal CDEFAULT_MOMENTUM, CDEFAULT_LEARNINGRATE;
	/** Create a multi-layer network.
	  * @param inputs The number of inputs taken in by the network
	  */
	#ifndef CONTROL
		MultiLayerNetwork(uint inputs);
	#else
		MultiLayerNetwork(uint inputs, uint neuronLabelOffset=0, PublicValues &pv=defaultControl);
		MultiLayerNetwork(uint inputs, PublicValues &pv);
	#endif

	/// Copy constructor, NOT YET IMPLEMENTED
	MultiLayerNetwork(MultiLayerNetwork &srcNet);

	/** Loads a multi-layer network from the given filename.
	  * @param filename The name of the file from which to load the network
	  * @throws Exception On any error
	  */
	MultiLayerNetwork(const std::string &filename);

	virtual ~MultiLayerNetwork();

	/// Adds a layer of the given size to the network. This new layer becomes the output layer
	virtual void addLayer(int size);

	/** Completely connects the given layer with the layer below it, i.e., all
	  * neurons in the given layer will give their output to all the neurons in the
	  * layer below. Weights are random.
	  * @param srcLayer The index of the source layer
	  * @throws Exception If the layer index provided is invalid
	  */
	virtual void connectLayer(uint srcLlayer);

	/** Connects two neurons, with random weight
	  * @param srcLayer The index of the layer in which the source neuron lies
	  * @param srcNrn The index of the source neuron in the source layer
	  * @param destNrn The index of the destination neuron (will be in the layer with index srcLayer+1)
	  * @throws Exception On any invalid argument values
	  */
	virtual void connect(uint srcLlayer, int srcNrn, int destNrn);

	/** Connects two neurons, with the weight provided.
	  * @param srcLayer The index of the layer in which the source neuron lies
	  * @param srcNrn The index of the source neuron in the source layer
	  * @param destNrn The index of the destination neuron (will be in the layer with index srcLayer+1)
	  * @param weight The weight of the link
	  * @throws Exception On any invalid argument values
	  */
	virtual void connect(uint srcLlayer, int srcNrn, int destNrn, real weight);

  virtual void connect(uint src, int srcNrn, uint dLayer, int destNrn, real weight);

	/** Sets the bias of the given neuron.
	  * @param layer The index of the layer in which the neuron lies
	  * @param nrn The index of the neuron in the layer
	  * @param bias The bias to be given to the neuron
	  */
	virtual void setBias(uint layer, int nrn, real bias);

	/// The number of layers in the network (does not count the input layer as a layer)
	virtual uint getLayerCount() const;

	/** Returns the output of the network for the given input.
	  * @param input A vector of getDimension() reals
	  * @return The corresponding output of the network
	  */
	virtual Vector getOutput(const Vector &input);

	/** Wrapper function to allow getOutput() to work for an array
	  * of real as input as well.
	  * Does exactly the same thing as Network::getOutput(real*).
	  */
	virtual Vector getOutput(real *input);

	/** Trains the network with data from the given TrainingSet using the
	  * backpropagation algorithm.
	  * @param T The TrainingSet containing input/desired-output vector pairs
	  * @param epochs The number of epochs to train the network. An epoch is a single
	  *				iteration through all input/desired-output vector pairs in T.
	  * @param learningRate The learning rate to be used for weight updation
	  * @param momentum The momentum factor to be used during weight updation. 0 by default.
	  */
	virtual void train(TrainingSet &T, uint epochs, real learningRate=DEFAULT_LEARNINGRATE, real momentum = DEFAULT_MOMENTUM);
	
	virtual void train(TrainingSet &T, Creal epochs, Creal learningRate=CDEFAULT_LEARNINGRATE, Creal momentum = CDEFAULT_MOMENTUM);
	
	/// get epochs, learningRate and momentum from the supplied PublicValues
	virtual void train(TrainingSet &T, PublicValues &parameters);

	/// compute error for the given traning set. The "epoch error" and "normalized epoch error" are then stored in the _control
	void getError(TrainingSet &ts);
	
	/// compute error for the given traning set. The "epoch error" and "normalized epoch error" are then stored in the _control
	void getErrorGREN(TrainingSet &ts);

//	struct Error;	
        struct  Error {
                /// [0] squared (the actual error function being minimalized)
                /// [1] absolute value of delta
                void operator+= (const Error &e)        { _val += e._val; }
                Error(const Vector &diff);
                Error() :_val(2) { zero(); }
                Error(real delta);

                /// publish current vals as name
                /// @param per if !0, also "normalized" version will be published
                void publish(std::string name, real per=0, PublicValues &pv=defaultControl);
                void zero() { _val[0] = 0; _val[1] = 0; }
                real getSq() const { return _val[0]; }
                real getAbs() const { return _val[1]; }
          protected:
                Vector _val;
        };
	/// real error for one example
	Error getErrorGREN(const Vector &input, const Vector &desired);

	/**
	 * Trains one example only
	 * @return error for this example
	 */
	virtual void trainExample(const Vector& input, const Vector& desiredOutput, real learningRate=DEFAULT_LEARNINGRATE, real momentum =  DEFAULT_MOMENTUM);

	/** Saves the network to the given filename.
	  * The file is a simple text file, open it up in a text editor
	  * to see the format, quite simple
	  * @param filename The name of the file to save the network in.
	  */
	virtual void save(const std::string &filename);

	/** Sets the activation function used by the neurons in the provided layer.
	  * @param layer The layer whose activation function is to be changed.
	  *				 layer>0 (as input neurons don't have any activation function)
	  *				 and layer<getLayerCount()
	  * @param f The activation function to be used
	  * @param df The derivative of the activation function to be used. Required for training.
	  * @throws Exception if an invalid layer is given
	  */
	virtual void setActivationFunction(uint layer, ActivationFunction f, ActivationFunction df);

	/// Returns "MultiLayerNetwork"
	virtual const char *getClassName() const;

	
	/// initialize weights to small values
	void resetWeights();

	/**
	 * Get a layer of the network (0=input)
	 * @throws Exception if an invalid layer is given
	 */
	const Layer &getLayer(uint layer) const;

	/**
	 * Get the total count of neurons in all layers (excluding input "layer")
	 */
	uint getNeuronsCount() const;

	/**
	 * Get count of all weights of all neurons
	 */
	uint getLinksCount() const;

	
	/**
	 * Warning: using this non-const version, you can change the network's behaviour.
	 * It's here because Neuron::getOutput() is not (yet) const
	 */
	Layer &getLayer(uint layer);

	/// get brief info about the topology, etc.
	operator std::string() const;
	#ifdef CONTROL
	PublicValues &getControl() { return *_control; }
	void setControl(PublicValues &ctrl) { _control = &ctrl; }

	void setLabelOffset(uint firstLabel) { _neuronLabelOffset = firstLabel; }

	/**
	 * Train the network using the error network.
	 * Generalized Relief Error Network - included in the MultiLayerNetwork, because it only provides new functions and can (better: must) be used together with MultiLayerNetwork
	 * The error network has getInputCount() + getOutputCount() inputs (in that order)
	 * @param ts contains input patterns to be trained. May also contain outputs - in this case, "real output error" will be computed
	 */
// smd -- removed as unneeded by ATHENA
	//virtual void trainGREN(MultiLayerNetwork &errorNetwork, TrainingSet &ts, uint epochs, real learningRate=DEFAULT_LEARNINGRATE, real momentum = DEFAULT_MOMENTUM);
  protected:
	/// attach or detach the gren. TODO: detach doesn't yet restore the input links of GREN, so it's unusable after training standalone
	void _attachGREN(MultiLayerNetwork &errorNetwork, bool detach);	

	///unlike _trainExample,  _trainExampleGREN assumes that input was already presented
	void _trainExampleGREN(MultiLayerNetwork &errorNetwork, const Vector& input, real learningRate, real momentum);

	PublicValues *_control;

	///Error summator
	///will be moved up for general use..
/*	struct	Error {
		/// [0] squared (the actual error function being minimalized)
		/// [1] absolute value of delta
		void operator+= (const Error &e)	{ _val += e._val; }
		Error(const Vector &diff);
		Error() :_val(2) { zero(); }
		Error(real delta);

		/// publish current vals as name
		/// @param per if !0, also "normalized" version will be published
		void publish(std::string name, real per=0, PublicValues &pv=defaultControl);
		void zero() { _val[0] = 0; _val[1] = 0; }
		real getSq() const { return _val[0]; }
		real getAbs() const { return _val[1]; }
	  protected:
		Vector _val;
	};
*/
	Error _exampleError; 
	#endif
	uint _neuronLabelOffset;	//start of labels of this network's neurons (--> 2 networks can interact)
};

///default MLP NumberParameters for parseArgs
#define DEFAULT_MLP_NPARS	\
	{ "learningRate", "learing rate (alpha)", MultiLayerNetwork::DEFAULT_LEARNINGRATE },	\
	{ "momentum", "learning momentum", MultiLayerNetwork::DEFAULT_MOMENTUM},	\
	{ "epochs", "how many epochs to train. You can usually keep it inf and interrupt by ESC if neccesary", 1000. },

}; //namespace annie
#endif // define _MULTILAYERNETWORK_H

