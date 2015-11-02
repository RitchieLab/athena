/**
 * Usable ValueUpdateListener implementations
 * @author OP
 */

#ifndef LISTENERS_H
#define LISTENERS_H

#include "Control.h"
//#include "cthread.h"

namespace annie	{

/**
 * Simple text-mode visualiser, which prints updated value after each change (synchronous).
 */
struct RigidVisualiser : ValueUpdateListener	{
	virtual void valueChanged(const Value &val);
};

/**
 * Simple text-mode visualiser, which prints updated value after each change (synchronous), it it didn't occur too soon
 */
struct SimpleVisualiser : RigidVisualiser {
	virtual void valueChanged(const Value &val);
};


//depends on cthread...
/*
struct Redrawer : ValueUpdateListener	{
	Redrawer(const char *triggeringValue) : _triggeringValue(triggeringValue) {}
	virtual void valueChanged(const Value &val)	{
		if(val.name() == _triggeringValue)
			forceRedraw();
	}
  protected:	
	const char *_triggeringValue;
};
*/

/**
 * Waits for keypress on each "epoch" if "stepped" is true
 */
struct Stepper: ValueUpdateListener	{
	virtual void valueChanged(const Value &val)	{
//		if(val.name() == "epoch" && !!control["stepped"])
//			waitForKey();
	}
};

class TrainingSet;
/**
 * Shuffles the given training set every nth epoch
 */
struct Shuffler	: public ValueUpdateListener{
	Shuffler(TrainingSet &ts, Creal shufflePeriod) : _ts(ts),  _shufflePeriod(shufflePeriod)	{}
	virtual void valueChanged(const Value &val);
  protected:
	TrainingSet &_ts;
	Creal _shufflePeriod;
};


}	//annie

#endif	//H
