#include "Exception.h"
#include <cmath>
#include "Vector.h"
#include "auxf.h"
#include <sstream>

using namespace std;
namespace annie
{
const Vector Vector::ZOID(0);	//0-dimensional constant

bool operator==(const Vector &o1, const Vector &o2) { 
	return ((vector<real> &)o1) == (vector<real> &) o2;
}
#define FORS	for(uint i=0; i<size(); i++)

std::ostream &operator<<(std::ostream &os, const annie::Vector &v)	{
	v.toStream(os);
	return os;
}

void Vector::toStream(std::ostream &s) const	{
	s << "(";
	if(size())	{
		uint i;
		for(i=0; i<size()-1; i++) s << (*this)[i] << ",";
		s << (*this)[i];
	}
	s << ")";
}

Vector Vector::fromInts(vector<int> iv)	{
	Vector o(iv.size());
	for(uint i=0; i<iv.size(); i++)
		o[i] = iv[i];
	return o;
}

Vector Vector::fromInts(uint size, int *data)	{
	Vector o(size);
	for(uint i=0; i<size; i++)
		o[i] = data[i];
	return o;
}

Vector::operator std::string() const	{
	stringstream ret;
	toStream(ret);
	return ret.str();
}
real Vector::magnitude()	const {
	real sum=0.;
	FORS	{
		real r = (*this)[i];
		sum += r * r;
	}
	return sqrt(sum);
}

void Vector::setAll(real val)	{
	FORS (*this)[i] = val;
}


void Vector::clamp(real min, real max)	{
	ASSERT(min <= max);
	FORS {
		real r = (*this)[i];
		if(r < min) (*this)[i] = min;
		else if(r > max) (*this)[i] = max;
	}
}

void Vector::normalize()	{
	(*this) *= (1./magnitude());
}

Vector Vector::operator- (const Vector &v) const	{
	Vector res(size());
	FORS	{
		res[i] = (*this)[i] - v[i];
	}
	return res;
}

Vector &Vector::operator+= (const Vector &v)	{
	FORS	{
		(*this)[i] += v[i];
	}
	return *this;
}

Vector &Vector::operator*= (const real r)	{
	FORS	{
		(*this)[i] *= r;
	}
	return *this;
}

Vector &Vector::operator+= (const real r)	{
	FORS	{
		(*this)[i] += r;
	}
	return *this;
}


Vector Vector::operator* (const real r)	const {
	Vector res = *this;
	res *= r;
	return res;
}


Vector Vector::operator+ (const real r) const	{
	Vector res = *this;
	res += r;
	return res;
}

Vector Vector::operator+ (const Vector &v)	const {
	Vector res = *this;
	res += v;
	return res;
}


real Vector::distance(const Vector &to) const	{
	Vector diff= (*this) - to;
	return diff.magnitude();
}

Vector &Vector::append(const Vector &v)	{
	for(uint i=0; i<v.size(); i++) push_back(v[i]);	//OPT: would InputIterator be faster?
	return *this;
}

void Vector::setMore(const Vector &v, uint startIndex)	{
	if(size() < startIndex + v.size()) resize(startIndex + v.size());
	for(uint i=0; i<v.size(); i++)
		(*this)[startIndex + i] = v[i];
}

Vector Vector::subset(uint start, uint size) const	{
	return Vector(begin() + start, begin() + start + size);
}

}; //namespace annie

