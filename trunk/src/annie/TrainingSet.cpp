#include "TrainingSet.h"
#include "Exception.h"
#include "File.h"
#include "Neuron.h"
#include "defines.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "examples.h"
#include "auxf.h"

using namespace std;
namespace annie
{

TrainingSet::TrainingSet(uint in,uint out)
{
	this->_nInputs=in;
	this->_nOutputs=out;
}

TrainingSet::~TrainingSet()
{}

void
TrainingSet::load_text(const string &filename)
{
	File file;
	try
	{
		file.open(filename);
	}
	catch (Exception &e)
	{
		string error(getClassName());
		error = error + "::" + getClassName() + "() - " + e.what();
		throw Exception(error);
	}

	string s;
	s=file.readWord();
	if (s.compare(getClassName()))
	{
		string error(getClassName());
		error = error + "::" + getClassName() + "() - File provided isn't a TrainingSet TEXT_FILE.";
		throw Exception(error);
	}
	while(!file.eof())
	{
		s=file.readWord();
		if (!s.compare("INPUTS"))
			_nInputs=file.readInt();
		else if (!s.compare("OUTPUTS"))
			_nOutputs=file.readInt();
		else if (!s.compare("IO_PAIRS"))
		{
			uint j;
			Vector input,output;

			while (!file.eof())
			{
				input.clear();
				output.clear();
				for (j=0;j<_nInputs;j++)
					input.push_back(file.readDouble());
				for (j=0;j<_nOutputs;j++)
					output.push_back(file.readDouble());
				_inputs.push_back(input);
				_outputs.push_back(output);
			}
		}
	}
}


void
TrainingSet::load_binary(const string &filename)
{
	ifstream file;
	double version;
	file.open(filename.c_str(),ios::binary);
	if (!file)
		throw Exception("TrainingSet::load_binary() - Couldn't open file for reading");
	file.read((char*)&version,sizeof(version));
	if (version!=atof(ANNIE_VERSION))
		throw Exception("TrainingSet::load_binary() - Invalid training set file encoutered (invalid version)");
	file.read((char*)&_nInputs,sizeof(_nInputs));
	file.read((char*)&_nOutputs,sizeof(_nOutputs));
	_inputs.clear();
	_outputs.clear();
	Vector v;
	real tmp;
	while (!file.eof())
	{
		v.clear();
		for (uint i=0;i<_nInputs;i++)
		{
			file.read((char*)&tmp,sizeof(tmp));
			v.push_back(tmp);
		}
		//Check this!! Why should it be giving EOF on read late?
		if (file.eof())
			break;
		_inputs.push_back(v);
		v.clear();
		for (uint i=0;i<_nOutputs;i++)
		{
			file.read((char*)&tmp,sizeof(tmp));
			v.push_back(tmp);
		}
		_outputs.push_back(v);
	}
	file.close();
}

TrainingSet::operator std::string() const	{
	return string(getClassName()) + ", " + int(getInputSize()) + "-->" + int(getOutputSize()) + ", " + int(_inputs.size()) + " pairs";
}

TrainingSet::TrainingSet(const string &filename, int file_type)
{
	_nInputs=_nOutputs==0;

	if (file_type == annie::TEXT_FILE)
		load_text(filename);
	else if (file_type == annie::BINARY_FILE)
		load_binary(filename);
	//else error
}

void
TrainingSet::addIOpair(real *input, real *output)
{
	Vector in,out;
	uint i;
	for (i=0;i<_nInputs;i++)
		in.push_back(input[i]);
	for (i=0;i<_nOutputs;i++)
		out.push_back(output[i]);
	addIOpair(in,out);
}

void
TrainingSet::addIOpair(const Vector &input, const Vector &output)
{
	ASSERT(input.size() == getInputSize());
	ASSERT(output.size() == getOutputSize());
	_inputs.push_back(input);
	_outputs.push_back(output);
}

void
TrainingSet::addIOpair(const Vector &input)
{
	addIOpair(input, Vector::ZOID);
}

bool
TrainingSet::epochOver() const
{
	if (_inputIter==_inputs.end() && _outputIter==_outputs.end())
		return true;
	return false;
}

void
TrainingSet::initialize()
{
	_inputIter=_inputs.begin();
	_outputIter=_outputs.begin();
}

void
TrainingSet::getNextPair(Vector &input, Vector &desired)
{
	if (_inputIter==_inputs.end())
	{
		string error(getClassName());
		error = error + "::getNextPair() - Passed the last I/O pair already. No more left.";
		throw Exception(error);
	}
	input=*_inputIter;
	desired=*_outputIter;
	_inputIter++;
	_outputIter++;
}

ostream& operator << (std::ostream& s, TrainingSet &T)
{
	Vector::iterator it;
	s<<T.getClassName()<<endl;
	if (s!=cout && s!=cerr)
		s<<"# TrainingSet information"<<endl;
	s<<"INPUTS "<<T._nInputs<<endl;
	s<<"OUTPUTS "<<T._nOutputs<<endl;
	if (s!=cout && s!=cerr)
		s<<"# -------------------------------------------------------- "<<endl;
	s<<"IO_PAIRS"<<endl;
	if (s!=cout && s!=cerr)
	{
		s<<"# -------------------------------------------------------- "<<endl;
		s<<"# Below follow lots of lines for each IO pair - a list of inputs"<<endl;
		s<<"# followed by a list of outputs"<<endl;
		s<<"# The first line will contain a vector with size INPUTS and"<<endl;
		s<<"# the next a vector of size OUTPUTS "<<endl;
	}
	vector < Vector >::iterator ioIn,ioOut;
	for (ioIn=T._inputs.begin(),ioOut=T._outputs.begin();ioIn!=T._inputs.end();ioIn++,ioOut++)
	{
		for (it=ioIn->begin();it!=ioIn->end();it++)
			s<<(*it)<<endl;
		s<<endl;
		for (it=ioOut->begin();it!=ioOut->end();it++)
			s<<(*it)<<endl;
		s<<endl;
		s<<endl;
	}
	return s;
}

void
TrainingSet::save_text(const string &filename)
{
	ofstream file;
	file.open(filename.c_str(),ios::out);
	if (!file)
		throw Exception("TrainingSet::save_text() - Couldn't open file for writing");
	file<<"ANNIE_FILE ";
	file<<ANNIE_VERSION;
	file<<endl;
	file<<"# Training Set information - the file integrity is"<<endl;
	file<<"# not checked when the file is loaded, so please do"<<endl;
	file<<"# not mess around with the file format as it may cause"<<endl;
	file<<"# errors that will be hard to trace"<<endl;
	file<<(*this);
	file.close();
}

void
TrainingSet::save_binary(const string &filename)
{
	ofstream file;
	file.open(filename.c_str(),ios::binary);
	if (!file)
		throw Exception("TrainingSet::save_binary() - Couldn't open file for writing");
	double version=atof(ANNIE_VERSION);
	file.write((char*)&version,sizeof(version));
	file.write((char*)&_nInputs,sizeof(_nInputs));
	file.write((char*)&_nOutputs,sizeof(_nOutputs));

	vector < Vector >::iterator ioIn,ioOut;
	Vector::iterator it;
	for (ioIn=_inputs.begin(),ioOut=_outputs.begin();ioIn!=_inputs.end();ioIn++,ioOut++)
	{
		for (it=ioIn->begin();it!=ioIn->end();it++)
			file.write((char*)&(*it),sizeof(*it));
		for (it=ioOut->begin();it!=ioOut->end();it++)
			file.write((char*)&(*it),sizeof(*it));
	}
	file.close();
}

void
TrainingSet::save(const string &filename, int file_type)
{
	if (file_type == TEXT_FILE)
		save_text(filename);
	else if (file_type == BINARY_FILE)
		save_binary(filename);
	else
	{
		string error(getClassName());
		error = error + "::save() - Invalid file type specified.";
		throw Exception(error);
	}
}


uint
TrainingSet::getSize() const
	{	return _inputs.size();	}

uint
TrainingSet::getInputSize() const
	{	return _nInputs;	}

uint
TrainingSet::getOutputSize() const
	{	return _nOutputs;	}

const char *
TrainingSet::getClassName() const
{	return "TrainingSet";	}

TrainingSet TrainingSet::operator+ (const TrainingSet &ts) const	{
	TrainingSet ret(*this);
	ret += ts;
	return ret;
}

TrainingSet &TrainingSet::operator+= (const TrainingSet &ts)	{
	if(getInputSize() != ts.getInputSize()) throw Exception("input sizes don't match!");
	if(getOutputSize() != ts.getOutputSize()) throw Exception("output sizes don't match!");

	TrainingSet b(ts);	//OP this 	is AWFUL - I copy the whole set only so that I can call initialize()... :)
	b.initialize();
	Vector in(getInputSize()), out(getOutputSize());
	while(!b.epochOver())	{
		b.getNextPair(in, out);
		addIOpair(in, out);
	}
	return (*this);
}

TrainingSet TrainingSet::xform(XformFunction ix, XformFunction ox, uint resI, uint resO)	{
	TrainingSet ret(resI, resO);
	Vector in(getInputSize()), out(getOutputSize());
	initialize();
	while(!epochOver())	{
		getNextPair(in, out);
		ret.addIOpair(ix(in), ox(out));
	}
	return ret;
}

TrainingSet TrainingSet::mixedXform(const TSTransformer &xf)	{
	uint resI = xf.getISize();
	uint resO = xf.getOSize();
	TrainingSet ret(resI, resO);
	Vector in1(getInputSize()), out1(getOutputSize());
	Vector in2(resI), out2(resO);
	initialize();
	while(!epochOver())	{
		getNextPair(in1, out1);
		xf.xform(in1, out1, in2, out2);
		ret.addIOpair(in2, out2);
	}
	return ret;
}

TrainingSet TrainingSet::xform(XformFunction ix, uint resI)	{
	return xform(ix, annie::vectorIdentity, resI, getOutputSize());
}

TrainingSet TrainingSet::xform(XformFunction ix)	{
	return xform(ix, getInputSize());
}


#if 0
void
TrainingSet::shuffle()
{
	int size = getSize()-1;
	vector< Vector >::iterator inIt,outIt;

	int chosen;
	while(size>=0)
	{
		chosen = (int)(fabs(random())*size);

		//inIt = &_inputs[chosen];
		//outIt = &_outputs[chosen];
		//OP: this seems portable
		inIt = _inputs.begin() + chosen;
		outIt = _outputs.begin() + chosen;

		_inputs.push_back(*inIt);
		_outputs.push_back(*outIt);

		inIt = _inputs.erase(inIt);
		outIt = _outputs.erase(outIt);
		size--;
	}
}

#else
/* OP: isn't that funny, I needed the shuffle() function so I made my own TrainingSet offspring. Later I decided to add it to the original TrainingSet and found this commented implementation in .cpp. :)
*
* my solution (uses some more memory but seems to be 3x faster - which is a bit surprising)
* I didn't put much effort measuring the difference. The shuffling is still very slow indeed. Looks like vector is too costly for large ammounts of data. I'd try to implement FieldTrainingSet, which would have arrays of in/out Vectors instead of vector<Vector>. The shuffling should be much faster then.
*/
void TrainingSet::shuffle()	{
	typedef vector<Vector> Vv;
	Vv ni, no;
	while(!_inputs.empty())	{
		assert(!_outputs.empty());
		int next = rand() % _inputs.size();
		Vv::iterator i = _inputs.begin() + next;
		Vector v = *i;
		_inputs.erase(i);
		ni.push_back(v);

		i = _outputs.begin() + next;
		v = *i;
		_outputs.erase(i);
		no.push_back(v);
	}

	_inputs.swap(ni);
	_outputs.swap(no);
}
#endif

}; //namespace annie

