#include "MultiLayerNetwork.h"
#include "Exception.h"
#include "SimpleNeuron.h"
#include "File.h"
#include "auxf.h"
#include "Control.h"
#include "math.h"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
namespace annie {
const real MultiLayerNetwork::DEFAULT_MOMENTUM = 0.3;
const real MultiLayerNetwork::DEFAULT_LEARNINGRATE = 0.2;
const Creal MultiLayerNetwork::CDEFAULT_LEARNINGRATE(defaultControl.init("defaultLearningRate", MultiLayerNetwork::DEFAULT_LEARNINGRATE));
const Creal MultiLayerNetwork::CDEFAULT_MOMENTUM(defaultControl.init("defaultMomentum", MultiLayerNetwork::DEFAULT_MOMENTUM));

//MultiLayerNetwork::Error::Error(real delta) :_val(2) { _val[0] = delta * delta; _val[1] = abs(delta ); }
// smd change to fabs from abs (check that it is ok)
MultiLayerNetwork::Error::Error(real delta) :_val(2) { _val[0] = delta * delta; _val[1] = fabs(delta ); }

MultiLayerNetwork::Error::Error(const Vector &diff) :_val(2) {
	zero();
	for(uint i=0; i<diff.size(); i++)
		(*this) += Error(diff[i]);
}

void MultiLayerNetwork::Error::publish(std::string name, real per, PublicValues &pv)	{
	pv.init(name, getAbs());
	if(per )	pv.init("normalized " + name, getAbs() / per);
}

#ifndef CONTROL
MultiLayerNetwork::MultiLayerNetwork(uint inputs) : Network(inputs,0), _neuronLabelOffset(0){
#else
MultiLayerNetwork::MultiLayerNetwork(uint inputs, uint neuronLabelOffset, PublicValues &pv) : Network(inputs,0), _control(&pv), _neuronLabelOffset( neuronLabelOffset) {
#endif
	_nLayers = 0;
	InputLayer *inputLayer = new InputLayer( _neuronLabelOffset,inputs);
	_layers.push_back(inputLayer);
}

//TODO:C&P
MultiLayerNetwork::MultiLayerNetwork(uint inputs, PublicValues &pv) : Network(inputs,0), _control(&pv), _neuronLabelOffset(0) {
	_nLayers = 0;
	InputLayer *inputLayer = new InputLayer( _neuronLabelOffset,inputs);
	_layers.push_back(inputLayer);
}


MultiLayerNetwork::MultiLayerNetwork(MultiLayerNetwork &src) : Network(src), _control(&defaultControl), _neuronLabelOffset(0) {
	throw Exception("MultiLayerNetwork::MultiLayerNetwork() - Copy constructor not yet implemented.");
}

MultiLayerNetwork::MultiLayerNetwork(const string &filename) : Network(0,0), _control(&defaultControl), 
		_neuronLabelOffset(0) {//TODO:read/save the idStart..
	_nLayers = 0;
	File file;
	try {
		file.open(filename.c_str());
	} catch (Exception &e) {
		string error(getClassName());
		error = error + "::" + getClassName() + "() - " + e.what();
		throw Exception(error);
	}

	string s;
	s=file.readWord();
	if (s.compare(getClassName())!=0) {
		string error(getClassName());
		error = error + "::" + getClassName() + "() - File supplied is not about this type of network.";
		throw Exception(error);
	}

	int maxLayerSize = Layer::MAX_LAYER_SIZE;
	while (!file.eof()) {
		s=file.readWord();
		if (!s.compare("INPUTS")) {
			_nInputs=file.readInt();
			InputLayer *inputLayer=new InputLayer(0,_nInputs);
			_layers.push_back(inputLayer);
		} else if (!s.compare("OUTPUTS"))
			_nOutputs=file.readInt();
		else if (!s.compare("MAX_LAYER_SIZE"))
			maxLayerSize=file.readInt();
		else if (!s.compare("SIZE"))
			addLayer(file.readInt());
		else if (!s.compare("Biases")) {
			vector < Layer * >::iterator l;
			for (l=_layers.begin()+1;l!=_layers.end();l++) {
				Layer *layer = (Layer *)(*l);
				for (uint j=0;j<layer->getSize();j++) {
					if (file.eof())
						break;
					SimpleNeuron &n = (SimpleNeuron&)layer->getNeuron(j);
					n.setBias(file.readDouble());
				}
			}
		} else if (!s.compare("BEGIN_META_DATA")) {
			static const basic_string <char>::size_type npos = (basic_string <char>::size_type)-1;
			string end("END_META_DATA");
			string metaData;
			s = file.readLine();
			while (s.find(end,0)==npos) {
				metaData = metaData + s + "\n";
				s = file.readLine();
			}
			if (metaData.length()>0)
				metaData.erase(metaData.length()-1); //getting rid of the last "\n"
			setMetaData(metaData);
		} else if (!s.compare("Connections"))
			break;
		else
			cerr<<getClassName()<<"::"<<getClassName()<<"() - Unrecognized token ("<<s<<"). Ignoring.\n";
	}
	//Now read connection information
	uint label1,label2;
	uint layer1,layer2;
	uint neuron1,neuron2;
	uint nLayers = getLayerCount();

	Layer *fromLayer,*toLayer;
	while (!file.eof()) {
		label1=file.readInt();
		file.readChar();
		label2=file.readInt();
		file.readChar();
		real wt=file.readDouble();

		layer1=label1/maxLayerSize;
		layer2=label2/maxLayerSize;

		neuron1=label1%maxLayerSize;
		neuron2=label2%maxLayerSize;

		if (layer1<0 || layer2<=0 || layer1>=nLayers || layer2>nLayers) {
			string error(getClassName());
			error = error + "::" + getClassName() + "() - Error in file. Label of neuron corresponds to an invalid layer.";
			throw Exception(error);
		}
		if (layer2-layer1!=1) {
			string error(getClassName());
			error = error + "::" + getClassName() + "() - Error in file. Connection not betwen ith and (i+1)th layer.";
			throw Exception(error);
		}

		fromLayer=_layers[layer1];
		toLayer=_layers[layer2];

		if (neuron1 >= fromLayer->getSize() || neuron2>=toLayer->getSize()) {
			string error(getClassName());
			error = error + "::" + getClassName() + "() - Error in file. There aren't that many neurons in the layer.";
			throw Exception(error);
		}
		connect(layer1,neuron1,neuron2,wt);
	} // while (!file.eof())
	file.close();
}

MultiLayerNetwork::~MultiLayerNetwork() {
	vector<Layer *>::iterator it;
	while (!_layers.empty()) {
		it = _layers.begin();
		Layer *l = (Layer *)(*it);
		delete l;
		it = _layers.erase(it);
	}
}

uint
MultiLayerNetwork::getLayerCount() const {
	return _nLayers;
}

Vector
MultiLayerNetwork::getOutput(real *input) {
	return Network::getOutput(input);
}

Vector
MultiLayerNetwork::getOutput(const Vector &input) {
	if (getLayerCount()==0) {
		string error(getClassName());
		error = error + "::getOutput() - There is no output layer";
		throw Exception(error);
	}
	_inputLayer()->setInput(input);
	try {
		return _outputLayer()->getOutput();
	} catch (Exception &e) {
		string error(getClassName());
		error = error + "::getOutput() - " + e.what();
		throw Exception(e);
	}
}

void
MultiLayerNetwork::addLayer(int size) {
	_nLayers++;
	int i;
	Layer *newLayer = new Layer(_neuronLabelOffset + _nLayers);
	for (i=0;i<size;i++)
		newLayer->addNeuron(new SimpleNeuron(newLayer->getLabel() * Layer::MAX_LAYER_SIZE + i));
	_layers.push_back(newLayer);
	_nOutputs = size;
}

void
MultiLayerNetwork::_connectLayer(Layer &srcLayer, Layer &destLayer)	{
	uint srcSize = srcLayer.getSize();
	uint destSize = destLayer.getSize();
	for (uint i=0;i<destSize;i++) {
		SimpleNeuron &destNeuron = (SimpleNeuron &)destLayer.getNeuron(i);
		for (uint j=0;j<srcSize;j++)
			destNeuron.connect(&srcLayer.getNeuron(j));
	}
}

void
MultiLayerNetwork::connectLayer(uint src) {
	if (src<0 || src>getLayerCount()-1) {
		string error(getClassName());
		error = error + "::connectLayer() - Invalid source layer specified.";
		throw Exception(error);
	}

	_connectLayer(*_layers[src], *_layers[src+1]);
}

void
MultiLayerNetwork::connect(uint src, int srcNrn, int destNrn, real weight) {
	Layer *srcLayer, *destLayer;
	uint srcSize, destSize;
	if (src<0 || src >= getLayerCount()) {
		string error(getClassName());
		error = error + "::connect() - Invalid source layer specified.";
		throw Exception(error);
	}

	srcLayer = _layers[src];
	destLayer = _layers[src+1];
	srcSize = srcLayer->getSize();
	destSize = destLayer->getSize();
	try {
		SimpleNeuron &destination = (SimpleNeuron&)destLayer->getNeuron(destNrn);
		destination.connect(&srcLayer->getNeuron(srcNrn),weight);
	} catch (Exception &e) {
		string error(getClassName());
		error = error + "::connect() - " + e.what();
		throw Exception(e);
	}
}

void
MultiLayerNetwork::connect(uint src, int srcNrn, uint dLayer, int destNrn, real weight) {
	Layer *srcLayer, *destLayer;
	uint srcSize, destSize;
	if (src<0 || src >= getLayerCount()) {
		string error(getClassName());
		error = error + "::connect() - Invalid source layer specified.";
		throw Exception(error);
  }
	if(dLayer > getLayerCount()){
		string error(getClassName());
		error = error + "::connect() - Invalid destination layer specified.";
		throw Exception(error);	  
	}

	srcLayer = _layers[src];
	destLayer = _layers[dLayer];
	srcSize = srcLayer->getSize();
	destSize = destLayer->getSize();
	try {
		SimpleNeuron &destination = (SimpleNeuron&)destLayer->getNeuron(destNrn);
		destination.connect(&srcLayer->getNeuron(srcNrn),weight);
	} catch (Exception &e) {
		string error(getClassName());
		error = error + "::connect() - " + e.what();
		throw Exception(e);
	}
}


void
MultiLayerNetwork::connect(uint srcLayer, int srcNrn, int destNrn) {
	connect(srcLayer,srcNrn,destNrn, AbstractNeuron::getRandomWeight());
}

void MultiLayerNetwork::resetWeights()	{
	for (uint i=1;i<=getLayerCount();i++) {
		uint layerSize = _layers[i]->getSize();
		for (uint j=0;j<layerSize;j++) {
			SimpleNeuron &n = (SimpleNeuron&)_layers[i]->getNeuron(j);
			n.randomizeWeights();
		}
	} //i = layers
}

void
MultiLayerNetwork::setBias(uint layer, int nrn, real bias) {
	if (layer<0 || layer>getLayerCount()) {
		string error(getClassName());
		error = error + "::setBias() - Invalid layer specified.";
		throw Exception(error);
	}
	try {
		SimpleNeuron &n = (SimpleNeuron&)_layers[layer]->getNeuron(nrn);
		n.setBias(bias);
	} catch (Exception &e) {
		string error(getClassName());
		error = error + "::setBias() - " + e.what();
		throw Exception(e);
	}
}

InputLayer *
MultiLayerNetwork::_inputLayer() {
	return (InputLayer*)_layers[0];
}

Layer *
MultiLayerNetwork::_outputLayer() {
	return _layers[getLayerCount()];
}

void MultiLayerNetwork::save(const string &filename)	{
	ofstream s;
	s.open(filename.c_str(),ios::out);
	if (!s) {
		string error(getClassName());
		error = error + "::save() - Could not open file for writing.";
		throw Exception(error);
	}

	s<<"ANNIE_FILE "<<ANNIE_VERSION<<endl;
	s<<"# Network information, the next line identifies the type of network"<<endl;
	s<<"# Constructing a network from a file which doesn't have the following identifier"<<endl;
	s<<"# should result in an error "<<endl;
	s<<"# DO NOT MAKE FORMAT CHANGES TO THIS FILE"<<endl;
	s<<"MultiLayerNetwork"<<endl;
	s<<"MAX_LAYER_SIZE "<<Layer::MAX_LAYER_SIZE<<" # There cannot be more than these many neurons in a layer"<<endl;
	s<<"INPUTS "<<getInputCount()<<endl;
	s<<"OUTPUTS "<<getOutputCount()<<endl;
	s<<"\nBEGIN_META_DATA"<<endl;
	if (!getMetaData().empty())
		s<<getMetaData()<<endl;
	s<<"END_META_DATA\n"<<endl;
	s<<"# Below follow layer sizes, the last one is the same as OUTPUTS"<<endl;
	s<<"# -------------------------------------------------------------------------"<<endl;

	vector < Layer * >::iterator it;
	uint ctr;
	for (ctr=0,it=_layers.begin()+1;it!=_layers.end();it++,ctr++)
		s<<"SIZE "<<(*it)->getSize()<<" # of layer "<<ctr<<endl;

	s<<"# Now follows information on the neural connections"<<endl;
	s<<"# Basically, MAX_LAYER_SIZE puts an upper bound on the number of neurons in a layer"<<endl;
	s<<"# so connections are specified as (label 1),(label 2),(weight) where "<<endl;
	s<<"# there is a connection from (label 1) --> (label 2) of weight (weight)"<<endl;
	s<<"# The label of neurons is MAX_LAYER_SIZE*L+N where L is the layer number and N is the"<<endl;
	s<<"# neuron number in that layer"<<endl;
	s<<"# Label of input layer is 0, layer below it is 1 and so on"<<endl;

	s<<"Biases"<<endl;
	for (ctr=0,it=_layers.begin()+1;it!=_layers.end();it++,ctr++) {
		s<<"# Layer "<<ctr<<", "<<(*it)->getSize()<<" lines follow"<<endl;
		for (uint i=0;i<(*it)->getSize();i++) {
			SimpleNeuron &n = (SimpleNeuron&)((*it)->getNeuron(i));
			s<<n.getBias()<<endl;
		}
	}

	s<<"Connections"<<endl;
	Neuron *nrn;
	for (it=_layers.begin();it!=_layers.end();it++) {
		uint layerSize = (*it)->getSize();
		for (ctr=0;ctr<layerSize;ctr++) {
			nrn=&((*it)->getNeuron(ctr));
			vector<int> labels;
			Vector weights;
			uint conns;
			conns=nrn->getInputs(labels,weights);
			for (uint j=0;j<conns;j++)
				s<<labels[j]<<","<<nrn->getLabel()<<","<<weights[j]<<endl;
		}
	}
	s.close();
}

void MultiLayerNetwork::trainExample(const Vector& input, const Vector& desired, real learningRate, real momentum)	{
	_inputLayer()->setInput(input);

	_exampleError.zero();	//error for this example... \sum (ex_i - o_i)^2
	uint nOutputs = getOutputCount();
	uint nLayers = getLayerCount();

	// set desired outputs
	for (uint i=0;i<nOutputs;i++) {
		SimpleNeuron &n = (SimpleNeuron&)_outputLayer()->getNeuron(i);
		n.setDesiredOutput(desired[i]);
		real delta = n.getOutput() - desired[i];
		_exampleError += Error(delta);
	}

	// back-propagate err - compute delta's
	for (int i=(int)nLayers;i>0;i--) {
		uint layerSize = _layers[i]->getSize();
		for (uint j=0; j<layerSize; j++) {
			SimpleNeuron &n = (SimpleNeuron&)_layers[i]->getNeuron(j);
			n.calculateNewWeights(learningRate,momentum);
		}
	} //i = layers

	// update weights
	for (uint i=1;i<=nLayers;i++) {
		uint layerSize = _layers[i]->getSize();
		for (uint j=0;j<layerSize;j++) {
			SimpleNeuron &n = (SimpleNeuron&)_layers[i]->getNeuron(j);
			n.update();
		}
	} //i = layers
}

void MultiLayerNetwork::train(TrainingSet &ts, Creal epochs, Creal learningRate, Creal momentum)	{
	if (ts.getInputSize() != getInputCount() || ts.getOutputSize() != getOutputCount()) {
		stringstream error;
		error << getClassName();
		error << "::train() - Training set not compatible with the network. (" << getInputCount() << "->" << getOutputCount() << " to be trained by : " << ts.getInputSize() << "->" << ts.getOutputSize();
		throw Exception(error.str());
	}
	ASSERT(epochs);

	try {
		Creal terminate = getControl().get("terminate");
		Creal epoch = getControl().init("epoch", 0);
		for (; epoch < epochs && !terminate; ++epoch) {
			Error epochErrorTmp;	//error summed for all examples in epoch
			ts.initialize();
			Vector input,desired;
			uint examples=0;
			while(!ts.epochOver()) {
				ts.getNextPair(input,desired);
				trainExample(input, desired, learningRate, momentum);
				epochErrorTmp += _exampleError;
				++examples;
			} // while ts
			 epochErrorTmp.publish("epoch error", examples, getControl());		} // epoch
		if(terminate ) cout << "terminated at epoch " << epoch;
	} catch (Exception &e) {
		string error(getClassName());
		error = error + "::train() - " + e.what();
		throw Exception(error);
	}
}

void
MultiLayerNetwork::train(TrainingSet &ts, uint epochs, real learningRate, real momentum) {
	PublicValues pv;	//private --> no listeners, etc.
	train(ts, pv.init("epochs", epochs), pv.init("learningRate", learningRate), pv.init("momentum", momentum));
}

void MultiLayerNetwork::getError(TrainingSet &t)	{
	train(t, 1, 0, 0);
}

void MultiLayerNetwork::train(TrainingSet &ts, PublicValues &p)	{
	train(ts, p["epochs"], p["learningRate"], p["momentum"]);
}

void
MultiLayerNetwork::_layerValid(uint layer) const {
	if (layer <= 0 || layer > _nLayers) {
		string error(getClassName());
		error = error + "::setActivationFunction() - Invalid layer provided";
		throw Exception(error);
	}
}

void
MultiLayerNetwork::setActivationFunction(uint layer, ActivationFunction f, ActivationFunction df) {
	_layerValid(layer);

	int N = _layers[layer]->getSize();
	for (int i=0; i<N; i++) {
		SimpleNeuron &n = (SimpleNeuron&)_layers[layer]->getNeuron(i);
		n.setActivationFunction(f,df);
	}
}

const char *
MultiLayerNetwork::getClassName() const {
	return "MultiLayerNetwork";
}


const Layer &MultiLayerNetwork::getLayer(uint layer) const {
	_layerValid(layer);
	return *_layers[layer];
}

Layer &MultiLayerNetwork::getLayer(uint layer) {
	_layerValid(layer);
	return *_layers[layer];
}

uint MultiLayerNetwork::getNeuronsCount() const	{
	uint c=0;
	for(uint i=1; i<=getLayerCount(); i++) c += getLayer(i).getSize();
	return c;
}

uint MultiLayerNetwork::getLinksCount() const	{
	uint c=0;
	for(uint i=1; i<=getLayerCount(); i++)	{
		const Layer &la = getLayer(i);
		for(uint n=0; n<la.getSize(); n++)	
			c += la.getNeuron(n).getInputCount();
	}
	return c;
}

MultiLayerNetwork::operator string() const	{
	stringstream res; res << getClassName() << " (" << getInputCount();
	for(uint i=1; i<=getLayerCount(); i++) res << " x " << getLayer(i).getSize();
	res << "), " << getLayerCount() << " layer(s), " << getNeuronsCount() << " neurons";
	return res.str();
}
}
; //namespace annie

