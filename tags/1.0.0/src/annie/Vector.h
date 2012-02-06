#ifndef Vector_H
#define Vector_H

#include "defines.h"
#include <vector>
//#include <iostream>
#include <iosfwd>

namespace annie
{
///A vector of real numbers
/* note: we should NOT use all of the std::vector functions. This opens us a way to define Vector as real[fixed] or dynamically resizable array of reals
 * TODO: hide them in the vector's descendant, then..!
 * TODO: it's inconvenient to set the size... but we cannot slow down [] ...
 * with all respect and understanding, I'd rename Vector to Vector
 * 
The time-critical function is [].
* OPT: employ STL's valarray ?
*/
//typedef std::vector<real> Vector;
typedef std::vector<real> VFather;
class Vector : protected  VFather	{
	public:
	#define BLUB_TYPE(t) typedef VFather::t t
	BLUB_TYPE(size_type);

	//TODO: purge dependency on iterators.... we can make iterated sucessor if needed. Or can we do it. without overhead ???
	BLUB_TYPE(iterator);
	BLUB_TYPE(const_iterator);
	///unhide selected functions...
	#undef BLUB_TYPE
	Vector()	{}
	///filled constructor
	Vector(uint items, real values)	: VFather(items, values) {}
	Vector(const Vector &src) : VFather(src)	{}
	Vector(size_type size) : VFather(size) {}
	Vector &operator = (const Vector &v)	{ VFather::operator=(v); return *this; }

	friend 	bool operator==(const Vector &o1, const Vector &o2);

	static Vector fromInts(std::vector<int> iv);
	static Vector fromInts(uint size, int *data);
	
	size_type size() const	{ return VFather::size(); }
	void reserve(size_type s) { VFather::reserve(s); }
	void resize(size_type s) { VFather::resize(s); }
	iterator begin() { return VFather::begin(); }
	const_iterator begin() const	{ return VFather::begin(); }
	const_iterator end() const	{ return VFather::end(); }
	real &operator[](uint i)	{ return VFather::operator[](i); }
	const real &operator[](uint i)	const { return VFather::operator[](i); }
	//unwanted?
	void clear() { VFather::clear(); }
	void push_back(const real & x)	{ VFather::push_back(x); }

	///more sugar here...
	real magnitude()	const;
	void normalize();
 	Vector operator- (const Vector &v) const;
	Vector operator+ (const Vector &v) const;
	Vector &operator+= (const Vector &v);

	Vector subset(uint start=0)	const {	return subset(start, size() - start); }
	Vector subset(uint start, uint size=0) const;

	///set all elements to val
	void setAll(real val);
	
	///set all elements to 0
	void null()	{ setAll(0.); }

	/**
	 * Insert contents of another vector at the end of this one.
	 * @return reference to self
	 */ 
	Vector &append(const Vector &v);	
	
	/**
	 * Will be resized if necessary
	 */
	void setMore(const Vector &v, uint startIndex=0);

	///all elements will be clamped to the given interval (inclusive)
	void clamp(real min, real max);
	
	///scalars
	Vector &operator*= (const real r);
	Vector operator* (const real r) const;
	Vector operator+ (const real r) const;
	Vector &operator+= (const real r);
	
	real distance(const Vector &to) const;
	const static Vector ZOID;	//0-dimensional constant
	
	/// example: "(1.4,4.4,4.5)"
	operator std::string() const;

	void toStream(std::ostream &s) const;
  protected:
	template<class InputIterator>
	Vector(InputIterator f,  InputIterator t) : VFather(f, t) {}
};
std::ostream &operator<<(std::ostream &os, const annie::Vector &v);
bool operator==(const Vector &o1, const Vector &o2);
}
#endif //_H
