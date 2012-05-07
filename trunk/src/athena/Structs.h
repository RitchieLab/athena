//Structs.h

#ifndef _STRUCTS_H
#define	_STRUCTS_H

#include <vector>
#include <set>

struct optSymbol{
  string symbol;
  bool noNT;
};


typedef vector<optSymbol> symbVector;

enum LogType{
  LogNone,
  LogSummary,
  LogDetailed,
};

#endif
