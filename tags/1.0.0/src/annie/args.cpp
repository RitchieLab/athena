/**
 * Dirty hack..
 */
#include "args.h"
#include "Control.h"
#include <map>
#include <stdlib.h>
#include <string.h>

using namespace std;
namespace annie	{

void parseArgs(int argc, char *argv[], const NumberParameter *numeric, StringParameter *string)	{
	ArgParser a(numeric, string);
	a.parse(argc, argv);
	a.print();
}

ArgParser::ArgParser(const NumberParameter *numeric, StringParameter *string, PublicValues &ctrl) : control(ctrl)	{
	uint i=0;
	for(const NumberParameter *n = numeric; n->name; n++, i++) {
		mn[n->name] = n;
		control.change(n->name, n->def);
	}

	i=0;
	for(StringParameter *n = string; n->name; n++, i++) {
		ms[n->name] = n;
		n->value = n->def;
	}
}

void ArgParser::parse(int argc, char *argv[])	{
	for(int j=1; j<argc; j++)	{
		char *a = argv[j];
		char *eq;
		if((eq = strchr(a, '=')))	{
			const char *val = eq+1;
			string name = string(a, 0, eq - a);
			if(ms.find(name) != ms.end()) ms[name]->value = val;
			else if(mn.find(name) == mn.end()) error(string("Dont know parameter '") + name + "'");
			else	{
				char *r;
				real rval = strtod(val, &r);
				if(r == a)	error(string("some mismatch near '") + a + "'");
				control.change(name, rval);
			}
		}	else	{	//=less
		    if(ms.find(a) != ms.end()) error("String param must have a value..");
		    else if(mn.find(a) != mn.end()) control.change(a, 1.);
			    else error(string("Dont know parameter '") + a + "'");
				}
			}
}

void ArgParser::print()	{
	for(Mn::const_iterator i=mn.begin(); i != mn.end(); i++)	{
		const NumberParameter *n = i->second;
		cout << "\t" << n->name << ": " << control[n->name] << endl;
	}

	for(Ms::const_iterator i=ms.begin(); i != ms.end(); i++)	{
		const StringParameter *n = i->second;
		cout << "\t" << n->name << ": '" << n->value << "'\n";
	}
}

void ArgParser::error(const string s)	{
	cout << "error: " << s << endl;
	help();
	throw Exception(s);
}

void ArgParser::help()	{
	cout << "syntax:"<<endl;
	cout << "numeric params:" << endl;

	for(Mn::const_iterator i=mn.begin(); i != mn.end(); i++)	{
		const NumberParameter *n = i->second;
		cout << "\t" << n->name << ": " << n->desc << ", default: " << n->def << endl;
	}

	cout << "string params:" << endl;
	for(Ms::const_iterator i=ms.begin(); i != ms.end(); i++)	{
		const StringParameter *n = i->second;
		cout << "\t" << n->name << ": " << n->desc << ", default: " << n->def << endl;
	}
}

}	//annie
