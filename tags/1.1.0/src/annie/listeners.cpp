#include "listeners.h"
#include "TrainingSet.h"
#include <iostream>

using namespace std;

namespace annie	{
		
void SimpleVisualiser::valueChanged(const Value &v)
{
	if(!v.isVisible() || !v.enoughTimePassed()) return;
	RigidVisualiser::valueChanged(v);
}

void RigidVisualiser::valueChanged(const Value &v)
{
	if(!v.isVisible()) return;
	real d = v - v.getLast();
	const char *dir = "==";
	if(d > 0) dir = "/\\";
	else if(d < 0) dir = "\\/";
	cout << v.name() << ": " << v << " " << dir << "\n";
}


void Shuffler::valueChanged(const Value &val) 	{
	if(val.name() != "epoch") return;
	if(!((int) val.get()  %  (int) _shufflePeriod )) 
			_ts.shuffle();
}
} //annie
