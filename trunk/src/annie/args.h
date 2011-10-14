/**
 * Command-line (and possibly config file) parameter handling
 * @author OP
 * $Id: args.h,v 1.1 2004/06/16 10:53:30 opx Exp $
 */

#ifndef EXAMPLES_H
#define EXAMPLES_H

#include "defines.h"
#include <map>

namespace annie
{

/*struct Parameter	{
	const char *name, * desc;
};*/

struct NumberParameter /*: Parameter*/	{
	const char *name, * desc;
	real def;	//default value
};

///also holds the destination
struct StringParameter /*: Parameter	*/{
	const char *name, * desc;
	const std::string def;
	std::string value;
};

class ArgParser	{
public:
	/**
	 * Parses commanline, stores numeric pars in control and string pars in StringParameter.value. 
	 * Also sets defaults if not specified.
	 * @param numeric, string - list of accepted arguments terminated by .name=NULL
	 */
	ArgParser(const NumberParameter *numeric, StringParameter *string, PublicValues &ctrl=defaultControl);
	void parse(int argc, char *argv[]);
	
	///print current values of the parameters
	void print();
protected:
	void error(const std::string txt);
	void help();
	typedef std::map<std::string, const NumberParameter*> Mn;
	typedef std::map<std::string, StringParameter*> Ms;
	Mn mn;
	Ms ms;
	PublicValues &control;
};

///shortcut - you don't have to use the ArgParser yourself
void parseArgs(int argc, char *argv[], const NumberParameter *numeric, StringParameter *strin);

}	//annie

#endif
