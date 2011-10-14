/**
 * Classes for sharing values (variables). This is useful e.g. for appropriately displaying current network error,
 * but can also be used to change learning parameters on-the-run.
 *
 * @file
 * @author OP
 */

#include "defines.h"
#ifdef CONTROL
#ifndef CONTROL_H
#define CONTROL_H

#include <map>
#include <list>
#include <string>

//struct ltstr { bool operator()(const char* s1, const char* s2) const { return strcmp(s1, s2) < 0;}};

struct SDL_mutex;
namespace annie
{
typedef std::string ValueKey;
//typedef const char *ValueKey;
class PublicValues;	//FW
class Value;	//FW
/**
 * Exported value. Used internally in PublicValues but Visualisers also access it.
 */
class Value	{
public:
	enum { VISIBLE=1, DEFAULT=VISIBLE};
	unsigned flags;	//<! we'd like to concentrate most of the flags used in listeners so that Listener doesn't have to hold per'value data and also that the flags are portable between listeners
	typedef real val_t;
	static const double DISPLAY_THRESHOLD;	//< value won't be re-displayer more often then this (seconds)

	void set(val_t to);
	val_t get() const { return _current; }
	val_t getLast() const { return _last; }
	Value &operator= (val_t val)	{ set(val); return *this; }
	val_t operator()() const	{ return get(); }
	operator val_t () const	{ return get(); }

	bool isVisible() const { return flags & VISIBLE; }
	Value(const std::string &name);

	/**
	 * It may be usable to track also processor time (clock(), etc.).
	 * That's why we use structure instead of type..
	 *
	 * TODO: copy TimeVal
	 */
	struct Time	{
		///sets all times to actual
		Time();
		/**
		 * Get difference from <code>t</code> in seconds.
		 * (The real-time component)
		 */
		double diffReal(const Time &t) { Time tmp = (*this) - t; return tmp.toSec(); }
		Time operator-(const Time &t) const;
		Time &operator-=(const Time &t)	{ (*this) = (*this) - t; return (*this); }

		double toSec() const;
		
		/**
		 * sets to current time
		 */
		void current();

		timeval _realCalendar;
	  protected:
		Time(const timeval &realCalendar) : _realCalendar(realCalendar) {}
	};
	std::string name() const	{ return _name; }
	
	/**
	 * Should be called when the visualiser doesn't want to produce too much (/too consuming) output.
	 * If true, the lastDisplayed time is reset to now
	 * TODO: only worx w/ one visualiser  (the visualisers themselves should hold this data)
	 */
	bool enoughTimePassed() const;

	double timeFromLastChange() const	{
		Time current;
		current -= _time;
		return current.toSec();
	}
	
private:
	real _current,	//< current value
	_last;		//< previous value
	Time _time, _lastTime;	//<time when current/last were set
	mutable Time	_lastTimeDisplayed;
	const std::string _name;
};

/**
 * Class used to store run-time values.
 *
 * Normal classes (code) update values.
 * Visualisers (GUI) present them to the users.
 *
 * The door is open also the other direction: Code reacts on changes done by GUI PublicValueslers.
 */
class PublicReal;
struct ValueUpdateListener;
class PublicValues
{
	friend class PublicReal;
public:
	PublicValues();
	~PublicValues();
	/**
	 * Change value. 
	 * If the value is changed for the first time, it will be registered w/ default _flags
	 * If the value is different from the previous, notification will be sent to the listeners
	 */
	void change(const ValueKey &name, real newVal);

	void inc(const ValueKey &name);

	/**
	 * Set flags for given value.
	 * @param flags bitwise combination of Value:: enum flags
	 */
	void setFlags(const ValueKey &name, unsigned flags);

	/**
	 * Lookup a PublicReal for the given name. Unlike get, this method throws if the value doesn't yet exist.
	 */
	//const Value &operator[](const std::string &name) const throw(Exception) { return (*this)[name]; }
	PublicReal operator[](const std::string &name) throw(Exception);

	/**
	 * Acquire a PublicReal for the given name. It will be created if it doesn't yet exist
	 */
	PublicReal get(const std::string &name);

	/**
	 * Same as get, but also set a value
	 * */
	PublicReal init(const std::string &name, real init=0.);
	
	///add change event listener - called automatically from <code>Listener</code>s
	void addListener(ValueUpdateListener *l);

	///remove change event listener
	void removeListener(ValueUpdateListener *l);

	///ehm... TODO, get rid of these
	void addListener(ValueUpdateListener &l) { addListener(&l);	}
	void removeListener(ValueUpdateListener &l)	{ removeListener(&l); }

	///print all contained vars
	void triggerAll() const;
private:
	///hook a direct reference to a value - must be unhooked
	void hookValue(const Value &v);
	void unHookValue(const Value &v)	{ /*TODO*/}

	//direct change - the Value must be present in this PublicValues
	void change(Value &v, real newVal);

	/**
	 * Spreads updates to the listeners
	 */
	void _updated(const Value &v) const;

	/**
	 * add value if it doesn't yet exist
	 * Must be called with locked valueLock !
	 */
	Value &_exist(const ValueKey &name);

	typedef std::map<std::string, Value> Values;
	Values _values;
	typedef std::list<ValueUpdateListener *> Listeners;
	Listeners _listeners;
#ifdef THREADS_ENABLED
	SDL_mutex *valuesLock;
#endif
};

PublicValues &getDefaultControl();
#define defaultControl	getDefaultControl()

/**
 * Callback interface for changing values
 *
 * @see listeners.h for some useful implementations
 */
struct ValueUpdateListener	{
	///registers listener to <code>control</code>
	ValueUpdateListener(PublicValues &ctrl=defaultControl);
	/**
	 * called whenever a value or it's _flags are changed
	 * Warning: must be reentrant!
	 */
	virtual void valueChanged(const Value &val) = 0;

	///helper if we move to multiple PublicValuess...
	PublicValues &getPublicValues() const;

	///deregisters listener
	virtual ~ValueUpdateListener();
  protected:
	PublicValues &control;
};



/**
 * This class has two purposes
 * - enable easy work w/ controlled variables
 * - enable transparent disabling of CONTROL framework (in that case, it behaves almost as an ordinary real &
 *
 * You can work with it almost as with "real &" (or real??), but keep in mind, that it can ANY TIME change the value under your hands
 */
//TODO: templatize parameter && enable other-than-real vars?
class PublicReal	{
public:
	friend class PublicValues;
	//TODO: also enable 'anon' (unnamed), so we can replace all reals eventually
	operator real()	const {	return get(); }
#ifdef FP_ENABLED
	int toInt() const { return 	get().toInt(); }
	double toDouble() const { return 	get().toDouble(); }
#endif
	//operator const real() const	{	return get(); }
	real operator=(real r) { set(r); return r; }
	real operator++ ()	{	set(get() + 1.); return get();	}
	real operator++(int)	{ real r = get(); set(r + 1); return r; }
	
	real operator-(real r) const { return get() - r; }
	real operator+(real r) const { return get() + r; }
	real operator*(real r) const { return get() * r; }
	
	real operator*=(real r) { set(get() * r); return get(); }
	real operator+=(real r) { set(get() + r); return get(); }
	real operator-=(real r) { set(get() - r); return get(); }

	double timeFromLastChange()	const { return _val.timeFromLastChange(); }

#ifdef CONTROL
	//we rely on never-lasting Values in control, which is not good. But we can count it on Value --> OK ..
	Value &_val;
	PublicValues &_control;
#else
	//ordinary real (+maybe name, if even..)
#endif
	~PublicReal()	{ _control.unHookValue(_val); }
	PublicReal(const PublicReal &p);
  private:
	void set(real r) { _control.change(_val, r); }
	real get() const { return _val.get(); }
	PublicReal(Value &v, PublicValues &c) : _val(v), _control(c) { _control.hookValue(_val); }
	PublicReal &operator=(const PublicReal &p);	//never implemented
};

typedef PublicReal Creal;

}; //namespace annie
#endif // define _EXCEPTION_H
#endif // CONTROL
