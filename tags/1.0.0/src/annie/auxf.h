/**
 * Varios auxiliary functions which are needed in almost every project ..
 * @author ??
 * $Id: auxf.h,v 1.1 2004/06/16 10:53:30 opx Exp $
 */
#ifndef AUXF_H
#define AUXF_H

#ifdef __cplusplus
extern "C"	{
#endif

#ifdef WIN32
	char *strdup(const char *s);
#else
#include <string.h>
#endif

	/**
	 * Versions throwing exceptions instead of NULL return
	 */
	void *xmalloc(int size);
	char *xstrdup(const char *s);
#ifdef __cplusplus

}	//export "C"

/**
 * backtrace support (+ demangle)
 * print backtrace to the logger (if supported on your platform)
 */
#include <string>
std::string debug_backtrace();
namespace annie	{
std::string operator+(const std::string &s, const int i);
inline std::string operator+(const std::string &s, const unsigned int i) { return s + (int) i; }
}

#endif	//C++

#endif	//_H
