/**
 * @author OP
 * $Id: Control.cpp,v 1.2 2004/06/16 10:48:34 opx Exp $
 */

#include "defines.h"
#ifdef CONTROL
#include "Control.h"
#include <sys/time.h>	//gettimeofday
#ifdef THREADS_ENABLED
#include <SDL/SDL_thread.h>
#endif

//#define const ValueKey &PublicValues::ValueKey
//XXX
using namespace std;
namespace annie
{
const double Value::DISPLAY_THRESHOLD=0.25;
PublicValues &getDefaultControl()	{
	static PublicValues ctrl;
	return ctrl;
}

ValueUpdateListener::ValueUpdateListener(PublicValues &ctrl) : control(ctrl)
{
//no more reasonable to register automatically!
//	getPublicValues().addListener(*this);
}

ValueUpdateListener::~ValueUpdateListener()
{
//	getPublicValues().removeListener(*this);
}

PublicValues &ValueUpdateListener::getPublicValues() const	{
	return defaultControl;
}

/*
   Taken from glibc-doc. Dunno why this is not in std library..

   Subtract the `struct timeval' _values X and Y,
   storing the result in RESULT.
   Return 1 if the difference is negative, otherwise 0.  */

int
timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y)
{
	/* Perform the carry for the later subtraction by updating y. */
	if (x->tv_usec < y->tv_usec) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_usec - y->tv_usec > 1000000) {
		int nsec = (x->tv_usec - y->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}

	/* Compute the time remaining to wait.
	   tv_usec is certainly positive. */
	result->tv_sec = x->tv_sec - y->tv_sec;
	result->tv_usec = x->tv_usec - y->tv_usec;

	/* Return 1 if result is negative. */
	return x->tv_sec < y->tv_sec;
}

inline double tv2sec(const timeval &tv)	{
	return tv.tv_sec + tv.tv_usec / (double) 1000000;
}

Value::Time::Time()	{
	current();
}

void Value::Time::current()	{
	gettimeofday(&_realCalendar, NULL);
}

Value::Time Value::Time::operator-(const Time &t) const	{
	timeval diff, b1=_realCalendar, b2=t._realCalendar;
	timeval_subtract(&diff, &b1, &b2);
	return Time(diff);
}

double Value::Time::toSec() const	{
	return tv2sec(_realCalendar);
}

bool Value::enoughTimePassed() const	{
	Time now;
	if(now.diffReal(_lastTimeDisplayed) >= DISPLAY_THRESHOLD)	{
		_lastTimeDisplayed = now;
		return true;
	} else return false;
}
// PublicValues
//////////////

PublicValues::PublicValues()
#ifdef THREADS_ENABLED
		:	valuesLock(SDL_CreateMutex())
#endif
{
}

PublicValues::~PublicValues()	{
#ifdef THREADS_ENABLED
	SDL_DestroyMutex(valuesLock);
#endif
}

void PublicValues::triggerAll() const	{
	for(Values::const_iterator i=_values.begin(); i!=_values.end(); i++)
		_updated(i->second);
}

PublicReal PublicValues::get(const string &name)	{
	#ifdef THREADS_ENABLED
	SDL_LockMutex(valuesLock);
	#endif
	Value &v = _exist(name);
	#ifdef THREADS_ENABLED
	SDL_UnlockMutex(valuesLock);
	#endif
	return PublicReal(v, *this);
}

PublicReal PublicValues::init(const string &name, real init)	{
	PublicReal r = get(name);
	r = init;
	return r;
}

void PublicValues::change(Value &v, real newVal)
{
	#ifdef THREADS_ENABLED
	SDL_LockMutex(valuesLock);
	#endif
	if(v.get() != newVal)	{
		v.set(newVal);
		#ifdef THREADS_ENABLED
		SDL_UnlockMutex(valuesLock);
		#endif
		_updated(v);
	} else	{
		#ifdef THREADS_ENABLED
		SDL_UnlockMutex(valuesLock);
		#endif
	}
}

void PublicValues::change(const ValueKey &name, real newVal)	{
	change(_exist(name), newVal);
}

void PublicValues::inc(const ValueKey &name)	{
	real v = (*this)[name].get();
	change(name, ++v);
}
/*
const Value	&PublicValues::operator[](const ValueKey &name) const throw(Exception) {
	#ifdef THREADS_ENABLED
	SDL_LockMutex(valuesLock);
	#endif
	Values::const_iterator i = _values.find(name);
	#ifdef THREADS_ENABLED
	SDL_UnlockMutex(valuesLock);
	#endif
	if(i == _values.end()) throw Exception(string("Value ") + name + " doesn't exist");
	return (*i).second;
}*/

PublicReal	PublicValues::operator[](const ValueKey &name) throw(Exception) {
	#ifdef THREADS_ENABLED
	SDL_LockMutex(valuesLock);
	#endif
	Values::iterator i = _values.find(name);
	#ifdef THREADS_ENABLED
	SDL_UnlockMutex(valuesLock);
	#endif
	if(i == _values.end()) throw Exception(string("Value ") + name + " doesn't exist");
	return PublicReal(i->second, *this);
}



void PublicValues::hookValue(const Value &v)	{
	//TODO
}

void PublicValues::setFlags(const ValueKey &name, unsigned flags)
{
	#ifdef THREADS_ENABLED
	SDL_LockMutex(valuesLock);
	#endif
	Value &v = _exist(name);
	v.flags = flags;
	#ifdef THREADS_ENABLED
	SDL_UnlockMutex(valuesLock);
	#endif

	_updated(v);
}

void PublicValues::_updated(const Value &v) const
{
	for(Listeners::const_iterator i=_listeners.begin(); i!=_listeners.end(); i++)
		(*i)->valueChanged(v);
}

Value &PublicValues::_exist(const ValueKey &name)	{
	Values::iterator i = _values.find(string(name));
	if(i == _values.end())	{
		_values.insert(Values::value_type(name, Value(name)));
		return _values.find(name)->second;
	} else return (*i).second;
}

void PublicValues::addListener(ValueUpdateListener *l)
{
	_listeners.push_front(l);
}

void PublicValues::removeListener(ValueUpdateListener *l)
{
	_listeners.remove(l);
}

Value::Value(const string &name) : flags(Value::DEFAULT), _current(0), _last(0), _lastTimeDisplayed(Time() - Time()), _name(name)	{ }

void Value::set(real to)	{
	_last = _current;
	_current = to;
	_lastTime = _time;
	_time.current();
}

/// PublicReal
/////////////

PublicReal::PublicReal(const PublicReal &p)	: _val(p._val), _control(p._control)	{
	_control.hookValue(_val);
}

}; //namespace annie
#endif // CONTROL
