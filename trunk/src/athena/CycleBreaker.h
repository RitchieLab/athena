/*
	Breaks cycles in a DAG using algorithm described in "A fast and effective heuristic for the feedback
arc set problem" Eades et. al.
*/

#ifndef __CYCLEBREAKER_H__
#define __CYCLEBREAKER_H__

#include <set>
#include <vector>
#include <map>
#include <deque>
#include <ga/GA2DBinStrGenome.h>

class CycleBreaker{

	public:

		int breakCycles(GA2DBinaryStringGenome& genome);


	private:
		struct node{
			std::set<int> parents;
			std::set<int> children;
		};

	bool removeSink(std::map<int,node>& network, std::deque<int>& s2);
	bool removeSources(std::map<int,node>& network, std::vector<int>& s1);
	void removeMaxSig(std::map<int,node>& network, std::vector<int>& s1);

};


#endif