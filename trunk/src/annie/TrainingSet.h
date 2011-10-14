#ifndef _TRAININGSET_H
#define _TRAININGSET_H

#include "Neuron.h"
#include "Network.h"

namespace annie
{

/// transforms one vector to another. Doesn't have to come up with the same size.
/// look in annie/examples.h for common xform functions
typedef Vector (*XformFunction)(const Vector &in);

/// 2 vectors -> 2 vectors
struct TSTransformer	{
	TSTransformer(uint isize, uint osize) : _isize(isize), _osize(osize)	{}
	virtual ~TSTransformer() {}

	/**
	 * @param out1, out2 are guaranted to be of get?Size() size. The function mustn't change the size.
	 */
	virtual void xform(const Vector &in1, const Vector &out1, Vector &in2, Vector &out2) const = 0;
	uint getISize() const { return _isize; }
	uint getOSize() const { return _osize; }
  private:
	uint _isize, _osize;
};

/** This is an abstraction for the set of patterns which are used to
  * "train" a network. It will contain sets of input and corresponding
  * desired-output vectors which will be used by the training algorithms
  * of supervised networks, such as the MultiLayerNetwork and RadialBasisNetwork
  *
  * You can save TrainingSets into a file. The file can be text, a simple format
  * so that you can construct the file from other programs as well, or the file
  * can be a binary file which is understood only by this library. Binary files
  * are generally smaller and quicker to load.
  *
  * You can use utilities provided with the library distribution to convert
  * from one format to the other.
  */
class TrainingSet
{
protected:
	/// The set of input vectors
	std::vector< Vector > _inputs;

	/// The set of corresponding desired output vectors
	std::vector< Vector > _outputs;

	/// An iterator through the input vectors
	std::vector< Vector >::iterator _inputIter;

	/// Iterator through the output vectors
	std::vector< Vector >::iterator _outputIter;

	/// Size of an input vector
	uint _nInputs;

	/// Size of an output vector
	uint _nOutputs;

	/// Save the training set in binary format to the given filename
	void save_binary(const std::string &filename);

	/// Save the file in text format to the given filename
	void save_text(const std::string &filename);

	/// Load from a binary file
	void load_binary(const std::string &filename);

	/// Load from a text file
	void load_text(const std::string &filename);
public:
	/** Create an empty training set
	  * @param in The size of an input vector
	  
	 * @param out The size of an output vector
	  */
	TrainingSet(uint in, uint out);

	/**
	 * Create a tranining set from a file.
	  * @param filename The name of the file
	  * @param file_type The type of the file (annie::BINARY_FILE) or (annie::TEXT_FILE). TEXT_FILE by default
	  * TODO: remake this constructor - TrainingSet(0,0) calls it :-)
	  */
	TrainingSet(const std::string &filename,int file_type = TEXT_FILE);

	virtual ~TrainingSet();

	/// Add an input and correspoding output vector
	void addIOpair(real *input,real *output);

	/** Add an input and corresponding output vector
	  * @throws Exception If the input vector size is not the same as getInputSize()
	  *					or the output vector size is not the same as getOutputSize()
	  */
	void addIOpair(const Vector &input, const Vector &output);

	///version without output - for unsipervised architectures.. Uses output dim 0
	void addIOpair(const Vector &input);

	/** Initialize the training set so that the first call to
	  * getNextPair() gives the first I/O pair stored.
	  * @see getNextPair
	  */
	virtual void initialize();

	/** Has a cycle through all input/output pairs been completed?
	  * @return true if all input/output pairs have been seen through getNextPair
	  * @see getNextPair
	  * @see initialize
	  */
	virtual bool epochOver() const;

	/// The number of input/desired-output vector pairs in the training set
	virtual uint getSize() const;

	/// The size of an input vector
	virtual uint getInputSize() const;

	/// The size of an output vector
	virtual uint getOutputSize() const;

	/** Returns the next input/output vector pair.
	  * You would typically use this in a fashion somewhat like:
	  * \code
	  * TrainingSet T("trset_file");
	  * T.initialize();
	  * Vector in,out;
	  * while (!T.epochOver())
	  * {
	  *		T.getNextPair(in,out);
	  *		// do what you need to with in and out
	  * }
	  * \endcode
	  * @param input The input vector is returned here
	  * @param desired The corresponding desired output vector is returned here
	  */
	virtual void getNextPair(Vector &input, Vector &desired);

	/// Allows you to print the TrainingSet onto a stream
	friend std::ostream& operator << (std::ostream& s, TrainingSet &T);

	/** Saves the training set to a file.
	  * @param filename The name of the file to save the training set to
	  * @param file_type The type of file to generate (annie::BINARY_FILE) or (annie::TEXT_FILE). Default is TEXT_FILE
	  */
	virtual void save(const std::string &filename, int file_type = TEXT_FILE);

	///debugging info - not to be confused with operator<< !
	operator std::string() const;
	
	/// Returns "TrainingSet"
	virtual const char *getClassName() const;

	/// Randomly changes order of stored pairs
	void shuffle();

	/**
	 * Weird TrainingSet arithmetics.
	 * Is not supposed to be effective: the only purpose is ease of use (at least for now)
	 * OPT: passing by value 
	 */
	TrainingSet operator+ (const TrainingSet &ts) const;
	TrainingSet &operator+= (const TrainingSet &ts);
	TrainingSet xform(XformFunction ix, XformFunction ox, uint resI, uint resO);
	TrainingSet xform(XformFunction ix, XformFunction ox)	{ return xform(ox, ix, getInputSize(), getOutputSize());	}

	///only the input part
	TrainingSet xform(XformFunction ix, uint resI);

	///keep dimensions
	TrainingSet xform(XformFunction ix);

	/**
	 * The transformer functions get 
	 */
	TrainingSet mixedXform(const TSTransformer &xf);
};

}; //namespace annie
#endif // define _TRAININGSET_H
