#include <string.h>
#include <stdlib.h>
#include "auxf.h"
#include <memory>	//bad_alloc
#include "defines.h"
#include <iostream>
#ifdef HAVE_BACKTRACE_SUPPORT
	#include <execinfo.h>
	#ifdef HAVE_CXXABI_H
	#include <cxxabi.h>
	#endif
#endif

#ifdef __cplusplus
using namespace std;

extern "C"	{
#define PISS_OFF	throw bad_alloc();
#else
	#define PISS_OFF	{ puts("not enough memory :-)"); abort();  }
#endif

#ifdef WIN32
	char *strdup(const char *s)	{
		int l = (int)strlen(s);

		char *r = (char *) malloc(l+1);
		if(!r) return NULL;
		memcpy(r, s, l+1);
		return r;
	}
#endif //WIN32

	void *xmalloc(int size)	{
		void *r = malloc(size);
		if(!r)	PISS_OFF
			return r;
	}

	char *xstrdup(const char *s)	{
		char *r = strdup(s);
		if(!r)	PISS_OFF
			return r;
	}



#ifdef __cplusplus
} //extern "C"

static string int2str(unsigned i)	{
	if(!i) return string("0");
	string s;
	while(i)	{
		s = (char)('0' + i%10) + s;	//slow :(
		i /= 10;
	}
	return s;
}

namespace annie	{
string operator+(const string &s, int i) {
	return s + int2str(i);
}
}

string debug_backtrace()	{
#ifdef HAVE_BACKTRACE_SUPPORT
	enum {ITEMS=50};
	void *array[ITEMS];
	uint size = backtrace (array, ITEMS);
	char **strings = backtrace_symbols (array, size);
	string tot = string("Obtained " ) + int2str(size) + " stack frames:\n";
	for (uint i = 0; i < size; i++)	{
#ifdef HAVE_CXXABI_H
		string all(strings[i]);
		uint s1 = all.find('(');
		uint lastslash=all.rfind('/', s1);
		uint s2 = all.find('+', s1);
		if(s1 == string::npos || s2 == string::npos ) tot += string("(No function name)  ") + all + "\n";
		else {
			string lib;
			if(lastslash != string::npos) lib =  all.substr(lastslash + 1, s1 - lastslash - 1);
			else lib = all.substr(0, s1+1);
			int status=-1;
			string mangled = string(all, s1+1, s2-s1-1);
			//log_assert.debug("from '" + all + "' extracted '" + mangled + "'");	//DD
			char *demangled = __cxxabiv1::__cxa_demangle(mangled.c_str(), NULL, 0, &status);
			string dm;
			//TODO: demangling a C-name results in bulshit (but status==OK). Is there a way to detect this?
			if(status)	tot += string("(DEMANGLE FAILED)  ") + mangled + " in " + lib + "\n";
			else tot += string(demangled) + " in " + lib + ' ' + all.substr(s2) + " (mangled: " + mangled + ")\n";
			if(demangled)	free(demangled);
		}
#else
		tot += string(strings[i]) + '\n';
#endif

	}
	if(size == ITEMS) tot += " (...) \n";
	free (strings);
	return tot;
#else
	return string("(no backtrace support enabled)");
#endif
}

#endif


