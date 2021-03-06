#include "vect3.h"

vect3::vect3(flt sx,flt sy,flt sz) : x (sx), y (sy), z(sz) {
	//cout << "Constructor" << endl;
}
vect3::vect3(const vect3 & a) : x (a.x), y (a.y), z(a.z) {
	//cout << "Constructor копирования" << endl;
}
vect3::~vect3(){
	//cout << "~Destructor" << endl;
}
flt vect3::operator*(const vect3 &a)const{
	return x * a.x + y * a.y + z * a.z;
}
vect3 vect3::operator*(flt l)const{
	return vect3(x*l, y*l, z*l);
}
vect3 vect3::operator/(flt l)const{
	return vect3(x / l, y / l, z / l);
}
bool vect3::operator==(const vect3 &a)const{
	return x == a.x && y == a.y && z == a.z;
}
bool vect3::operator!=(const vect3 &a)const{
	return x != a.x || y != a.y || z != a.z;
}
const flt &vect3::operator[](int i)const{
	if(i == 0)
		return x;
	if(i == 1)
		return y;
	if(i == 2)
		return z;
	throw 1;
}
flt &vect3::operator[](int i){
	if(i == 0)
		return x;
	if(i == 1)
		return y;
	if(i == 2)
		return z;
	throw 1;
}
vect3& vect3::operator=(const vect3& a){
	x = a.x;
	y = a.y;
	z = a.z;
	return *this;
}
vect3 &vect3::operator*=(flt l){
	x *= l;
	y *= l;
	z *= l;
	return *this;
}
vect3 &vect3::operator/=(flt l){
	x /= l;
	y /= l;
	z /= l;
	return *this;
}
vect3 vect3::operator-()const{
	return vect3(-x , -y, -z);
}
vect3& vect3::operator+=(const vect3 & a){
	x += a.x;
	y += a.y;
	z += a.z;
	return *this;
}
vect3& vect3::operator-=(const vect3 & a){
	x -= a.x;
	y -= a.y;
	z -= a.z;
	return *this;
}
void vect3::show()const{
	cout << x << "; " << y << "; " << z << endl;
}
vect3 vect3::operator+(const vect3 &a)const{
	return vect3(x+a.x , y+a.y, z+a.z);
}
vect3 vect3::operator-(const vect3 &a)const{
	return vect3(x-a.x , y-a.y, z-a.z);
}
vect3 operator*(flt l,vect3 const & a){
	return vect3(a.x * l, a.y * l, a.z * l);
}
ostream& operator<<(ostream& os,const vect3& a){
    os << a.x << "\t" << a.y << "\t" << a.z;
    return os;
}
istream& operator>>(istream& is, vect3& a){
    flt x,y,z;
    is >> x >> y >> z;
    return is;
}
flt vect3::mod() const{
	flt a = sqr();
	if (a == 0)
		return 0;
	return sqrt(a);
}
flt vect3::sqr() const{
	return *this * *this;
}

flt vect3::normMax() const{
	flt a = fabs(x), b = fabs(y), c = fabs(z);
	if(a > b)
		if(a > c)
			return a;
		else
			return c;
	else
		if(b > c)
			return b;
		else
			return c;
}

vect3 vect3::unit() const{
	flt m = mod();
	if(m == 0)
		return vect3();
	return *this / m;
}

vect3 vProduct(vect3 const &a, vect3 const &b) {
	return vect3(a.y*b.z - a.z*b.y, -a.x*b.z + a.z*b.x, a.x*b.y - a.y*b.x);
}

vect3 rotate(vect3 const &r, flt a, flt b, flt c) { //euler angles
	vect3 rs;
	rs.x = r * vect3(cos(a)*cos(c)-sin(a)*cos(b)*sin(c), -cos(a)*sin(c)-sin(a)*cos(b)*cos(c), sin(a)*sin(b));
	rs.y = r * vect3(sin(a)*cos(c)+cos(a)*cos(b)*sin(c), -sin(a)*sin(c)+cos(a)*cos(b)*cos(c), -cos(a)*sin(b));
	rs.z = r * vect3(sin(b)*sin(c), sin(b)*cos(c), cos(b));
	return rs;
}
