#ifndef _vect3_H_
#define _vect3_H_

#include<iostream>
#include<cmath>

typedef float flt;

using namespace std;

struct vect3{
	vect3(flt sx=0,flt sy=0, flt sz=0);
	vect3(const vect3 & a);
	~vect3();
	vect3 operator-()const;
	vect3 &operator+=(const vect3 & a);
	vect3 operator+(const vect3 &a)const;
	vect3 &operator-=(const vect3 & a);
	vect3 operator-(const vect3 &a)const;
	flt operator*(const vect3 &a)const;
	vect3 operator*(flt l)const;
	vect3 operator/(flt l)const;
	friend vect3 operator*(flt l,vect3 const & a);
	bool operator==(const vect3 &a)const;
	bool operator!=(const vect3 &a)const;
	const flt &operator[](int i)const;
	flt &operator[](int i);
	vect3 &operator=(const vect3& a);
	vect3 &operator*=(flt l);
	vect3 &operator/=(flt l);
	void show()const;
	friend ostream& operator<<(ostream& os,const vect3& a);
	friend istream& operator>>(istream& is, vect3& a);
	flt mod() const;
	flt sqr() const;
	flt normMax() const;
	vect3 unit() const;
	
	flt x, y, z;
};

vect3 vProduct(vect3 const &a, vect3 const &b);
vect3 rotate(vect3 const &r, flt a, flt b, flt c);

#endif  // _vect3_H_
