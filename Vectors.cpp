#include <iostream>
#include "Vectors.h"
#include <cstring>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>

using namespace std;

Point::Point()  // default constructor
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
	Point(x, y, z);
}

Point::Point(float x, float y, float z) : x(x), y(y), z(z) {}  // constructor

Point:: ~Point() {} // destructor

Point::Point(const Point& rhs) // copy constructor
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
}

Point& Point:: operator=(const Point& rhs) // operator=
{
	if (this != &rhs)
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
	}
	return *this;
}

bool Point:: operator==(Point& rhs) // if two points coincide
{
	return(x == rhs.x && y == rhs.y && z == rhs.z);
}

ostream& operator<<(ostream& out, const Point& rhs) // cout
{
	out << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
	return out;
}

//##########################################################################################################################

Vector::Vector(): Point() // default constructor
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
}

Vector::Vector(float x, float y, float z):Point(x, y, z), x(x), y(y), z(z){}

Vector::Vector(float x, float y, float z, Point* p1, Point* p2):Point(x, y, z) // constructor
{
	x = p2->x - p1->x;
	y = p2->y - p1->y;
	z = p2->z - p1->z;
}

Vector:: ~Vector() {} // destructor

Vector::Vector(const Vector& rhs): Point(rhs)// copy constructor
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
}

Vector& Vector:: operator=(const Vector& rhs) // operator=
{
	if (this != &rhs)
	{
	    Point::operator=(rhs);
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
	}
	return *this;
}

float Vector::vector_length()// finds length of vector
{
	return sqrt(x * x + y * y + z * z);
}

Vector Vector::vector_direction()// finds vector direction
{
	try {
		if (vector_length() == 0) {
			throw obj; // throws exception
		}
		else {
			Vector Single(x / vector_length(), y / vector_length(), z / vector_length());
			return Single;
		}
	}
	catch (VectorLengthException& p) //check if vector length is invalid, enter new vector
	{
		cerr << "Caught " << p.what() << endl;
		cout << "Please, enter another vector!" << endl;
		float x, y, z;
		cout << "\tEnter values for x, y, z to define the Vector V1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Vector V1(x, y, z);
	}
}

bool Vector::zero_vector_check() // checks if vector is zero
{
	return vector_length() == 0;
}

bool Vector::parallel_vectors(Vector& v) // checks if vectors are parallel
{
	try
	{
		if (vector_length() == 0 || v.vector_length() == 0)
		{
			throw obj; // throws exception if length equals 0
		}
		else if ((this->x / v.x) == (this->y / v.y) && (this->x / v.x) == (this->z / v.z) && (this->y / v.y) == (this->z / v.z))
		{
			return true;
		}
		else
			return false;

	}
	catch (VectorLengthException& p) { //checks if vector is invalid and if so user must enter another vector
		cerr << p.what() << endl;
		cout << "Please, enter another vector!" << endl;
		float x, y, z;
		cout << "\tEnter values for x, y, z to define the Vector V1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Vector V1(x, y, z);
	}
}

bool Vector::perpendicular_vectors(Vector& v) // check if vectors are perpendicular
{
	try
	{
		if (vector_length() == 0 || v.vector_length() == 0)
		{
			throw obj; // throws exception
		}

		if ((this->x * v.x) + (this->y * v.y) + (this->z * v.z) == 0)
		{
			return true;
		}

		else
			return false;

	}
	catch (VectorLengthException& p) {
		cerr << p.what() << endl;
		cerr << "Caught " << p.what() << endl;
		cout << "Please, enter another vector!" << endl;
		float x, y, z;
		cout << "\tEnter values for x, y, z to define the Vector V1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Vector V1(x, y, z);

	}

}

Vector Vector:: operator+(const Vector& v) // sum of two vectors
{
	Vector temp(x + v.x, y + v.y, z + v.z);
	return temp;
}

Vector Vector:: operator-(const Vector& v) // subtraction of two vectors
{
	Vector temp(x - v.x, y - v.y, z - v.z);
	return temp;
}

float Vector:: operator*(const Vector& v) // vector multiply of two vectors
{
	float temp = x * v.x + y * v.y + z * v.z;
	return temp;
}


Vector Vector:: operator*(int num) // multiply of vector and number
{
	Vector temp(x * num, y * num, z * num);
	return temp;
}

Vector Vector:: operator^(const Vector& v) // scalar multiply of two vectors
{
	Vector temp(y * v.z - z * v.y, -(x)*v.z + z * v.x, x * v.y - y * v.x);
	return temp;
}

float Vector:: operator()(Vector& v, Vector& w) // mixed multiply
{
	float vecSum;
	vecSum = (v.x * w.y - w.x * v.y) * this->z - (this->x * w.y - w.x * this->y) * v.z + (this->x * v.y - v.x * this->y) * w.z;
	return vecSum;
}

ostream& operator<<(ostream& out, const Vector& rhs)
{
	out << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
	return out;
}


//############################################################################################################################

Line::Line(): Vector() // default constructor
{
	starting_point.x = 0;
	starting_point.y = 0;
	starting_point.z = 0;
	line_direction = Vector(0, 0, 0);
}

Line::Line(float x, float y, float z, Point& p1, Vector& v): Vector(x, y, z) // constructor
{
	starting_point.x = p1.x;
	starting_point.y = p1.y;
	starting_point.z = p1.z;
	line_direction = v;
}

Line::Line(float x, float y, float z, Point& p1, Point& p2): Vector(x, y, z) // constructor
{
	Vector temp(p1, p2);
	starting_point.x = p1.x;
	starting_point.y = p1.y;
	starting_point.z = p1.z;
	line_direction = temp;
}

Line:: ~Line() {} // destructor

Line::Line(const Line& rhs):Vector(rhs) // copy constructor
{
	starting_point.x = rhs.starting_point.x;
	starting_point.y = rhs.starting_point.y;
	starting_point.z = rhs.starting_point.z;
	line_direction = rhs.line_direction;
}

Line& Line:: operator=(const Line& rhs) // operator=
{
	if (this != &rhs)
	{
	    Vector::operator=(rhs);
		starting_point.x = rhs.starting_point.x;
		starting_point.y = rhs.starting_point.y;
		starting_point.z = rhs.starting_point.z;
		line_direction = rhs.line_direction;
	}
	return *this;
}

Vector Line::find_direction() // finds direction of line
{
	return line_direction;
}

Vector Line::normal_vector() // finds the normal vector
{
	Vector normal;
	Vector V(starting_point.x, starting_point.y, starting_point.z);
	normal = line_direction ^ V;
	return normal;
}

float Line::Angel_between_two_lines(Line& rhs) // finds the angle between two lines
{
	float angle = (this->line_direction.x * rhs.line_direction.x +
		this->line_direction.y * rhs.line_direction.y +
		this->line_direction.z * rhs.line_direction.z) /
		(this->line_direction.vector_length() * rhs.line_direction.vector_length());
	return acos(angle);
}

bool Line:: operator+(Point& p) // checks if Point& p lays on the line
{
	return ((p.x - starting_point.x / line_direction.x) == (p.y - starting_point.y / line_direction.y) == (p.z - starting_point.z / line_direction.z));
}

bool Line:: operator||(Line& rhs) // checks if two lines are parallel
{
	Vector temp = this->line_direction ^ rhs.line_direction;
	return (temp.vector_length() == 0);
}

bool Line:: operator==(Line& rhs) // checks if implicit line coincides with Line& rhs
{
	return (this->starting_point.x == rhs.starting_point.x && this->starting_point.y == rhs.starting_point.y
		&& this->starting_point.z == rhs.starting_point.z && this->line_direction == rhs.line_direction);
}

bool Line:: operator|(Line& rhs) // checks if two lines are perpendicular
{
	return (this->line_direction.perpendicular_vectors(rhs.line_direction));
}


//#########################################################################################################################


Segment::Segment(): Line() // default constructor
{
	this->start_point = Point(0, 0, 0);
	this->end_point = Point(0, 0, 0);
}

Segment::Segment(Point& p1, Vector& v, Point& p2): Line(p1, v,) // constructor
{
	this->start_point = p1;
	this->end_point = p2;
	Line l(this->start_point, this->end_point);
}

Segment:: ~Segment() {} // destructor

Segment::Segment(const Segment& rhs):Line(rhs) // copy constructor
{
	this->start_point = rhs.start_point;
	this->end_point = rhs.end_point;
	Line l(this->start_point, this->end_point);
}

Segment& Segment:: operator=(const Segment& rhs) // operator=
{
	if (this != &rhs)
	{
	    Line::operator=(rhs);
		this->start_point = rhs.start_point;
		this->end_point = rhs.end_point;
		Line l(this->start_point, this->end_point);
	}
	return *this;
}

float Segment::find_length_of_Segment() // finds the length of the segment
{
	Vector temp;
	temp.x = this->end_point.x - this->start_point.x;
	temp.y = this->end_point.y - this->start_point.y;
	temp.z = this->end_point.z - this->start_point.z;
	return temp.vector_length();
}

Point Segment::find_Mid_point_of_Segment() // finds the middle point of the segment
{
	Point Min_point;
	Min_point.x = (this->start_point.x + this->end_point.x) / 2;
	Min_point.y = (this->start_point.y + this->end_point.y) / 2;
	Min_point.z = (this->start_point.z + this->end_point.z) / 2;
	return Min_point;
}

bool Segment:: operator==(Point& p) // checks if a point lays on the segment
{
	float x0 = this->start_point.x;
	float y0 = this->start_point.y;
	float z0 = this->start_point.z;

	float x1 = this->end_point.x;
	float y1 = this->end_point.y;
	float z1 = this->end_point.z;

	float x = p.x;
	float y = p.y;
	float z = p.z;

	float AB = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0));
	float AP = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0));
	float PB = sqrt((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y) + (z1 - z) * (z1 - z));

	return AB == AP + PB;
}


//#########################################################################################################################


Triangle::Triangle(): Point() // default constructor
{
	this->A = Point(0, 0, 0);
	this->B = Point(0, 0, 0);
	this->C = Point(0, 0, 0);
}

Triangle::Triangle(float x, float y, float z, Point& p1, Point& p2, Point& p3): Point(x, y, z) // constructor
{
	try
	{
		if (p1.x == p2.x == p3.x && p1.y == p2.y == p3.y && p1.z == p2.z == p3.z)
		{
			throw obj1; // throw exception
		}
		else
		{
			this->A = p1;
			this->B = p2;
			this->C = p3;
		}
	}
	catch (EqualPointException& p) {
		cerr << p.what() << endl;
		cerr << p1 << " " << p2 << " " << p3 << endl;
		cout << "\tEnter values for x, y, z to define the starting Point P1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P1(x, y, z);

		cout << "\tEnter values for x, y, z to define the ending Point P2" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P2(x, y, z);

		cout << "\tEnter values for x, y, z to define the ending Point P3" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P3(x, y, z);

		Triangle T(P1, P2, P3);


	}
}

Triangle:: ~Triangle() {} // destructor

Triangle::Triangle(const Triangle& rhs):Point(rhs) // copy constructor
{
	this->A = rhs.A;
	this->B = rhs.B;
	this->C = rhs.C;
}

Triangle& Triangle:: operator=(const Triangle& rhs) // operator=
{
	if (this != &rhs)
	{
	    Point::operator=(rhs);
		this->A = rhs.A;
		this->B = rhs.B;
		this->C = rhs.C;
	}
	return *this;
}

void Triangle::define_Type() // finds the type of the triangle(equilateral or acute)
{
	Vector AB(A, B);
	Vector AC(A, C);
	Vector BC(B, C);
	if (AB.vector_length() == AC.vector_length() == BC.vector_length())
		cout << " Triangle is Equilateral and is Acute " << endl; // Равностранен и остроъгълен
	else if ((AB.vector_length() == AC.vector_length() && AB.vector_length() != BC.vector_length()) ||
		(AB.vector_length() == BC.vector_length() && AB.vector_length() != AC.vector_length()) ||
		(AC.vector_length() == BC.vector_length() && AC.vector_length() != AB.vector_length()))
		cout << "Triangle is Isosceles and is Acute" << endl; // Равнобедрен и остроъгълен
	else if (AB.vector_length() != AC.vector_length() != BC.vector_length())
		cout << "Triangle is Scalene" << endl; // разностранен
}

double Triangle::find_Area() // finds area of the triangle
{
	Vector AB(A, B);
	Vector AC(A, C);
	Vector BC(B, C);
	double p = (AB.vector_length() + AC.vector_length() + BC.vector_length()) / 2;
	double area = sqrt(p * (p - AB.vector_length()) * (p - BC.vector_length()) * (p - AC.vector_length()));
	return area;
}

double Triangle::find_Perimeter() // find perimeter
{
	Vector AB(A, B);
	Vector AC(A, C);
	Vector BC(B, C);
	double p = (AB.vector_length() + AC.vector_length() + BC.vector_length());
	return p;
}

Point Triangle::find_Medicenter() // find medicenter
{
	Point m((A.x + B.x) / 2, (A.y + B.y) / 2, (A.z + B.z) / 2);
	Vector Am(C, m);
	Point medicenter((C.x + m.x) / 3, (C.y + m.y) / 3, (C.z + m.z) / 3);
	return medicenter;
}

float Triangle::does_point_lays_on_(Point& p) // checks if a point lays on the triangle
{
	Vector AB(A, B);
	Vector AC(A, C);

	AB^ AC;
	float D = (AB ^ AC).x * A.x + (AB ^ AC).y * A.y + (AB ^ AC).z * A.z; // A + Bx + Cx + ((D))
	return (AB ^ AC).x * p.x + (AB ^ AC).y * p.y + (AB ^ AC).z * p.z + D;
}

bool Triangle:: operator<(Point& p) // checks if Point& p lays on the surface of the triangle and is inner
{
	float min_x = (((A.x < B.x ? A.x : B.x) < C.x) ? (A.x < B.x ? A.x : B.x) : C.x);
	float max_x = (((A.x > B.x ? A.x : B.x) > C.x) ? (A.x > B.x ? A.x : B.x) : C.x);

	float min_y = (((A.y < B.y ? A.y : B.y) < C.y) ? (A.y < B.y ? A.y : B.y) : C.y);
	float max_y = (((A.y > B.y ? A.y : B.y) > C.y) ? (A.y > B.y ? A.y : B.y) : C.y);

	float min_z = (((A.z < B.z ? A.z : B.z) < C.z) ? (A.z < B.z ? A.z : B.z) : C.z);
	float max_z = (((A.z > B.z ? A.z : B.z) > C.z) ? (A.z > B.z ? A.z : B.z) : C.z);
	return (p.x > min_x && p.x < max_x && p.y > min_y && p.y < max_y && p.z > min_z && p.z < max_z);
}

bool Triangle:: operator>(Point& p) // checks if Point& p lays on the surface of the triangle and is outer
{
	float min_x = (((A.x < B.x ? A.x : B.x) < C.x) ? (A.x < B.x ? A.x : B.x) : C.x);
	float max_x = (((A.x > B.x ? A.x : B.x) > C.x) ? (A.x > B.x ? A.x : B.x) : C.x);

	float min_y = (((A.y < B.y ? A.y : B.y) < C.y) ? (A.y < B.y ? A.y : B.y) : C.y);
	float max_y = (((A.y > B.y ? A.y : B.y) > C.y) ? (A.y > B.y ? A.y : B.y) : C.y);

	float min_z = (((A.z < B.z ? A.z : B.z) < C.z) ? (A.z < B.z ? A.z : B.z) : C.z);
	float max_z = (((A.z > B.z ? A.z : B.z) > C.z) ? (A.z > B.z ? A.z : B.z) : C.z);

	if (p.x < min_x || p.x > max_x)
		if (p.y < min_y || p.y > max_y)
			if (p.z < min_z || p.z > max_z)
				return true;
}

bool Triangle:: operator==(Point& p) // if Point p lays on at least one line of the triangle
{
	Segment AB(A, B);
	Segment AC(A, C);
	Segment BC(B, C);

	return (AB.operator==(p)) || (AC.operator==(p)) || (BC.operator==(p));
}


//##########################################################################################################################


Tetrahedron::Tetrahedron() // default constructor
{
	this->p1 = Point();
	this->p2 = Point();
	this->p3 = Point();
	this->p4 = Point();
	Tetrahedron(p1, p2, p3, p4);
}

Tetrahedron::Tetrahedron(Point& p1, Point& p2, Point& p3, Point& p4) // checks if points coincide and constructs
{
	try {
		if (((p1.x == p2.x && p1.y == p2.y && p1.z == p2.z) ||
			(p1.x == p3.x && p1.y == p3.y && p1.z == p3.z) ||
			(p1.x == p4.x && p1.y == p4.y && p1.z == p4.z) ||
			(p2.x == p3.x && p2.y == p3.y && p2.z == p3.z) ||
			(p2.x == p4.x && p2.y == p4.y && p2.z == p4.z) ||
			(p3.x == p4.x && p3.y == p4.y && p3.z == p4.z)))
			throw obj1; // throws exception

		else {
			this->p1 = p1;
			this->p2 = p2;
			this->p3 = p3;
			this->p4 = p4;
		}
	}
	catch (EqualPointException& p) {

		cerr << p.whatt() << endl;
		cout << "Please initialize the Tetrahedron again!" << endl; //enter new values
		float x, y, z;
		int choice;
		cout << "\tEnter values for x, y, z to define the starting Point P1" << endl; //defining point one
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P1(x, y, z);

		cout << "\tEnter values for x, y, z to define the ending Point P2" << endl; //defining point two
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P2(x, y, z);

		cout << "\tEnter values for x, y, z to define the ending Point P3" << endl; //defining point three
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P3(x, y, z);

		cout << "\tEnter values for x, y, z to define the ending Point P4" << endl; //defining point four
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P4(x, y, z);

		Tetrahedron Te(P1, P2, P3, P4); //constructing the tetrahedron
	}
}


Tetrahedron:: ~Tetrahedron() {} // destructor

Tetrahedron::Tetrahedron(const Tetrahedron& rhs):Point(rhs) // copy constructor
{
	p1 = rhs.p1;
	p2 = rhs.p2;
	p3 = rhs.p3;
	p4 = rhs.p4;
}

Tetrahedron& Tetrahedron:: operator=(const Tetrahedron& rhs) // operator=
{
	if (this != &rhs)
	{
	    Point::operator=(rhs);
		p1 = rhs.p1;
		p2 = rhs.p2;
		p3 = rhs.p3;
		p4 = rhs.p4;
	}
	return *this;
}

bool Tetrahedron::check_if_Regular() // checks if all lines are equal
{
	Triangle a(p1, p2, p3);
	Triangle b(p1, p4, p3);
	Triangle c(p2, p4, p3);
	Triangle d(p1, p2, p4);
	if (a.find_Area() == b.find_Area() && b.find_Area() == c.find_Area() && c.find_Area() == d.find_Area()) //checks if all areas are equal
		return true;
	else
		return false;
}

bool Tetrahedron::check_if_Ortogonal() //checks if the tetrahedron is ortogonal
{
	Segment a(p1, p2);
	Segment b(p1, p3);
	Segment c(p2, p3);
	Segment d(p1, p4);
	Segment e(p2, p4);
	Segment f(p3, p4);

	return (a * a + f * f == b * b + e * e == c * c + d * d);
}

float Tetrahedron::get_Surface() // lateral area
{
	Triangle a(p1, p2, p3);
	Triangle b(p1, p4, p3);
	Triangle c(p2, p4, p3);
	Triangle d(p1, p2, p4);
	float res = a.find_Area() + b.find_Area() + c.find_Area() + d.find_Area();
	return res;
}

float Tetrahedron::get_Volume() // finds volume
{
	Vector a(p1, p2);
	Vector b(p1, p3);
	Vector c(p1, p4);
	float temp = a.operator()(b, c);
	return (temp / 6);
}

bool Tetrahedron:: operator<(Point& p) // checks if Point& p lays on the surface of the tetrahedron(on the surface of at least one of its triangles) and is inner
{
	Triangle ABC(p1, p2, p3);
	Triangle ABD(p1, p2, p4);
	Triangle BCD(p2, p3, p4);
	Triangle ACD(p1, p3, p4);
	return (ABC.operator<(p) || ABD.operator<(p) || BCD.operator<(p) || ACD.operator<(p));
}

bool Tetrahedron:: operator>(Point& p) // checks if Point& p lays on the surface of the tetrahedron(on the surface of at least one of its triangles) and is outer
{
	Triangle ABC(p1, p2, p3);
	Triangle ABD(p1, p2, p4);
	Triangle BCD(p2, p3, p4);
	Triangle ACD(p1, p3, p4);
	return (ABC.operator>(p) || ABD.operator>(p) || BCD.operator>(p) || ACD.operator>(p));
}

bool Tetrahedron:: operator==(Point& p) // checks if Point& p lays at least on one line of the tetrahedron
{
	Segment AB(p1, p2);
	Segment AC(p1, p3);
	Segment BC(p2, p3);
	Segment AD(p1, p4);
	Segment BD(p2, p4);
	Segment CD(p3, p4);

	return (AB.operator==(p)) || (AC.operator==(p)) || (BC.operator==(p) || AD.operator==(p)) || (BD.operator==(p)) || (CD.operator==(p));
}


//#######################################################################################################################


bool Menu::readData(char choose) // here we choose to read from console or to read from a file
{
	if (choose == '1')
		return true;
	else
		return false;
}


ofstream out("output.txt", ios::out);

void Menu::mainMenu() //Main menu function
{
	while (other_object)
	{
		cout << "   ********************************************************   " << endl;
		cout << "************************** M e n u  **************************" << endl;
		cout << "   ********************************************************   " << endl;
		cout << endl << endl << endl;
		cout << "1. Point" << endl;
		cout << "2. Vector" << endl;
		cout << "3. Line" << endl;
		cout << "4. Segment" << endl;
		cout << "5. Triangle" << endl;
		cout << "6. Tetrahedron" << endl << endl;

		other_operation = true;

		cout << "\tPress 1 if you want to enter from the console" << endl << endl;
		cout << "\tPress a number different than 1, if you want to enter from a file" << endl << endl;
		char ch;
		cin >> ch; //Based on this input we either read from the console or froam a file
		choose = ch;

		cout << endl;
		cout << "Please choose an option from above or press 0 to exit program" << endl;
		int choice;

		cin >> choice; //Reads user input about which operation should be done


		if (choice == 0) //Exits the program
		{
			break;
		}

		while (other_operation)
		{
			switch (choice) {
			case 1: //Loads point menu
				pointMenu();
				break;
			case 2: //Loads vector menu
				vectorMenu();
				break;
			case 3: //Loads line menu
				lineMenu();
				break;
			case 4: //Loads segment menu
				segmentMenu();
				break;
			case 5: //Loads triangle menu
				triangleMenu();
				break;
			case 6: //Loads vector menu
				tetrahedronMenu();
				break;
			}
		}
    }
}

void Menu:: pointMenu() //Point menu
{
	if (readData(choose)) //If user choosed to enter from the console
	{
	    out.open("output.txt");
		float x, y, z;
		int choice;
		cout << "\tEnter values for x, y, z to define the Point P1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P1(x, y, z);
		cout << "P1: " << P1 << endl << endl;
		out << "P1: " << P1 << endl << endl;

		while (other_operation) //Point menu operations
		{
			cout << "************************** P O I N T **************************" << endl;
			cout << endl << endl << endl;
			cout << "PLEASE, CHOOSE OPERATION FROM BELOW OR PRESS 0 TO GO BACK TO MAIN MENU" << endl << endl;
			cout << "\t1. Check if two points coincide" << endl << endl;
			cin >> choice; //Enter which operation should the program do

			if (choice == 0) //Goes back to main menu if input is 0
			{
				break;
			}

			switch (choice) { //Checks which operation has the user choosed and does the following action
			case 1: //Checks if two points coincide
			{
				cout << "\tEnter values for x, y, z to define Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P2(x, y, z);
				cout << "P2" << P2 << endl << endl;
				out << "P2" << P2 << endl << endl;

				if (P1.operator==(P2)){
					cout << "Points: P1" << P1 << " = P2" << P2 << endl << endl;
                    out << "Points: P1" << P1 << " = P2" << P2 << endl << endl;
				}
				else{
					cout << "Points: P1" << P1 << " != P2" << P2 << endl << endl;
					out << "Points: P1" << P1 << " != P2" << P2 << endl << endl;
				}
				break;
			}
			default:
            {
                //If input is invalid
				cout << "Invalid choice!" << endl;
                out << "Invalid choice!" << endl;
				break;
            }
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //While loop continues
			else {
				other_operation = false; //While loop breaks
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y') {
					other_object = true; //Returns the user to the main menu
				}
				else
					other_object = false; //Exits program
			}
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}
	else //Reading from a file
	{
		ifstream input("Point.txt", ios::in);//Opens file Point.txt
		if (!input) //If file is invalid throws error
		{
			cerr << "File cannot be opened!" << endl;
		}
		input.seekg(1, ios::beg);

		out.open("output.txt");

		float x, y, z;
		input >> x;
		input >> y;
		input >> z;

		Point P1(x, y, z);
		cout << "P1: " << P1 << endl;  //Prints point one values
		out << "P1: " << P1 << endl;

		int choice;
		input >> choice; //Reads which operation should be done

		switch (choice) {
		case 1: //Checks if two points coincide and prints the result
		{
			input >> x;
			input >> y;
			input >> z;

			Point P2(x, y, z);
			cout << "P2: " << P2 << endl << endl;
			out << "P2: " << P2 << endl << endl;
			input.close();

			if (P1.operator==(P2)){
				cout << "Points: P1: " << P1 << " = P2: " << P2 << endl << endl;
				out << "Points: P1: " << P1 << " = P2: " << P2 << endl << endl;
			}
			else{
				cout << "Points: P1: " << P1 << " != P2: " << P2 << endl << endl;
				out << "Points: P1: " << P1 << " != P2: " << P2 << endl << endl;
			}
			break;
		}
		default:
        {
		     //If input is invalid
			cout << "Invalid choice!" << endl;
			out << "Invalid choice!" << endl;
			break;
        }
    }
		other_operation = false;
		cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
		char c;
		cin >> c;
		if (c == 'Y') {
			other_object = true; //Returns user to the main menu
		}
		else
			other_object = false; //Exits program

		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}
}




void Menu::vectorMenu() //Vector menu
{
	if (readData(choose)) //If user choose to enter from the console
	{
	    out.open("output.txt");
		int choice;
		float x, y, z;
		cout << "\tEnter values for x, y, z to define the Vector V1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Vector V1(x, y, z);
		cout << "V1" << V1 << endl << endl;
		out << "V1" << V1 << endl << endl;

		bool flag = true;

		while (other_operation) //Vector menu operations
		{
			cout << "************************** V E C T O R **************************" << endl;
			cout << endl << endl << endl;
			cout << "PLEASE, CHOOSE OPERATION FROM BELOW OR PRESS 0 TO GO BACK TO MAIN MENU" << endl << endl;
			cout << "\t1.Find Vector length" << endl << endl;
			cout << "\t2.Find Vector direction" << endl << endl;
			cout << "\t3.Check if NULL Vector" << endl << endl;
			cout << "\t4.Check if two Vectors are parallel" << endl << endl;
			cout << "\t5.Check if two Vectors are perpendicular" << endl << endl;
			cout << "\t6.Sum of two Vectors" << endl << endl;
			cout << "\t7.Multiply Vector with a number" << endl << endl;
			cout << "\t8.Multiply Vector with a Vector" << endl << endl;
			cout << "\t9.Scalar multiplication of two Vectors" << endl << endl;
			cout << "\t10.Mixed multiplication of three Vectors" << endl << endl;

			cin >> choice;

			if (choice == 0) //Checks if user wants to go back to main menu and has entered 0
			{
				break;
			}

			switch (choice) { //Cheks which operation the user chose and does the following action
			case 1: // Calculates and prints vector length
			{
				cout << "\tVector length is: " << V1.vector_length() << endl << endl;
				out << "\tVector length is: " << V1.vector_length() << endl << endl;
				break;
			}
			case 2: //Finds and prints vector's direction
			{
				cout << "\tVector direction is: " << V1.vector_direction() << endl << endl;
				out << "\tVector direction is: " << V1.vector_direction() << endl << endl;
				break;
			}
			case 3: //Checks if vector is null and prints the result
			{
				if (V1.zero_vector_check())
                {
                    cout << "Vector " << V1 << " is NULL vector " << endl;
					out << "Vector " << V1 << " is NULL vector " << endl;
                }
				else
                {
                    cout << "Vector " << V1 << " IS NOT NULL vector " << endl;
                    out << "Vector " << V1 << " IS NOT NULL vector " << endl;
                }
				break;
			}
			case 4: //Checks if two vectors are parallel
			{
				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates second vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;
				if (V1.parallel_vectors(V2))// Prints result based on entered vector points
				{
					cout << "\tVectors " << V1 << " and " << V2 << " are parallel" << endl << endl;
					out << "\tVectors " << V1 << " and " << V2 << " are parallel" << endl << endl;
				}
				else
                {
 					cout << "\tVectors " << V1 << " and " << V2 << " are NOT parallel" << endl << endl;
 					out << "\tVectors " << V1 << " and " << V2 << " are NOT parallel" << endl << endl;
                }
				break;
			}
			case 5: //Checks if two vectors are perpendicular to one another
			{
				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates second vector
                cout << "Vector: " << V2 << "initialized!" << endl;
                out << "Vector: " << V2 << "initialized!" << endl;
				if (V1.perpendicular_vectors(V2)) //Prints result based on entered vector points
				{
					cout << "\tVectors " << V1 << " and " << V2 << " are perpendicular" << endl << endl;
					out << "\tVectors " << V1 << " and " << V2 << " are perpendicular" << endl << endl;

				}
				else
                {
					cout << "\tVectors " << V1 << " and " << V2 << " are NOT perpendicular" << endl << endl;
					out << "\tVectors " << V1 << " and " << V2 << " are NOT perpendicular" << endl << endl;
                }
				break;
			}
			case 6: //Sums two vectors
			{
				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates second vector
				cout << "Vector: " << V2 << "initialized!" << endl;
                out << "Vector: " << V2 << "initialized!" << endl;

				Vector V3 = V1 + V2; //Sums the two vectors
				cout << "Vector " << V1 << " + Vector " << V2 << " is equal to Vector " << V3 << endl << endl; //Prints third vector with the sum of

				out << "Vector " << V1 << " + Vector " << V2 << " is equal to Vector " << V3 << endl << endl;																								   //the two entered from the user
				break;
			}
			case 7: //Multiplies vector with a number
			{
				float x;
				cout << "\tx = " << endl;
				cin >> x;
				Vector V2 = V1 * x;
				cout << "Vector " << V1 << " multiplied with " << x << " is equal to Vector " << V2 << endl << endl; //Prints the result
				out << "Vector " << V1 << " multiplied with " << x << " is equal to Vector " << V2 << endl << endl; //Prints the result
				break;
			}
			case 8: //Multiplies vector with a vector
			{
				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates second vector
				cout << "Vector: " << V2 << "initialized!" << endl;
                out << "Vector: " << V2 << "initialized!" << endl;

                float result = V1 * V2;
				cout << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to " << result << endl << endl; //Prints the result
				out << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to " << result << endl << endl; //Prints the result
				break;
			}
			case 9: //Does scalar multiplication of two vectors
			{
				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates second vector
				cout << "Vector: " << V2 << "initialized!" << endl;
                out << "Vector: " << V2 << "initialized!" << endl;
                Vector V3 = V1 ^ V2;
				cout << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to Vector " << V3 << endl << endl; //Prints the result
				out << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to Vector " << V3 << endl << endl; //Prints the result
				break;
			}
			case 10: //Does mixed multiplication of vectors
			{
				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates second vector
				cout << "Vector: " << V2 << "initialized!" << endl;
                out << "Vector: " << V2 << "initialized!" << endl;

				cout << "\tEnter values for x, y, z to define Vector V3" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V3(x, y, z); //Creates third vector
				cout << "Vector: " << V3 << "initialized!" << endl;
                out << "Vector: " << V3 << "initialized!" << endl;

                float result = V1.operator()(V2, V3);
				cout << "Mixed multiplication of Vector " << V1 << " Vector " << V2 << " and Vector " << V3 << " is equal to " << result << endl << endl; // Prints the result
				out << "Mixed multiplication of Vector " << V1 << " Vector " << V2 << " and Vector " << V3 << " is equal to " << result << endl << endl; // Prints the result
				break;
			}
			default:
            {
			    cout << "Invalid choice!" << endl;
                out << "Invalid choice!" << endl;
				break;
            } //If input is invalid
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //While loop continues
			else {
				other_operation = false; //While loop breaks
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y')
					other_object = true; //Returns to main menu
				else
					other_object = false; //Exits program
			}
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}
	else
	{
		ifstream input("Vector.txt"); //Opens file vector.txt
		if (!input) //If file is invalid throws error
		{
			cerr << "File cannot opened!" << endl;
		}
		input.seekg(1, ios::beg);

		out.open("output.txt");
		float x, y, z;
		input >> x;
		input >> y;
		input >> z;

		Vector V1(x, y, z);
		cout << "V1: " << V1 << endl << endl; //Prints vector 1 points
		int choice;
		input >> choice;

		switch (choice) { //Calculates vector's length
		case 1: // Calculates and prints vector length
			{
				cout << "\tVector length is: " << V1.vector_length() << endl << endl;
				out << "\tVector length is: " << V1.vector_length() << endl << endl;
				break;
			}
			case 2: //Finds and prints vector's direction
			{
				cout << "\tVector direction is: " << V1.vector_direction() << endl << endl;
				out << "\tVector direction is: " << V1.vector_direction() << endl << endl;
				break;
			}
			case 3: //Checks if vector is null and prints the result
			{
				if (V1.zero_vector_check())
                {
                    cout << "Vector " << V1 << " is NULL vector " << endl;
					out << "Vector " << V1 << " is NULL vector " << endl;
                }
				else
                {
                    cout << "Vector " << V1 << " IS NOT NULL vector " << endl;
                    out << "Vector " << V1 << " IS NOT NULL vector " << endl;
                }
				break;
			}
		case 4: //Checks if two vectors are parallel
		{
			input >> x;
			input >> y;
			input >> z;
			Vector V2(x, y, z); //Creates second vector based on file information
			cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;
				if (V1.parallel_vectors(V2))// Prints result based on entered vector points
				{
					cout << "\tVectors " << V1 << " and " << V2 << " are parallel" << endl << endl;
					out << "\tVectors " << V1 << " and " << V2 << " are parallel" << endl << endl;
				}
				else
                {
 					cout << "\tVectors " << V1 << " and " << V2 << " are NOT parallel" << endl << endl;
 					out << "\tVectors " << V1 << " and " << V2 << " are NOT parallel" << endl << endl;
                }
				break;
		}
		case 5: //Checks if tow vectors are perpendicular
		{
			input >> x;
			input >> y;
			input >> z;
			Vector V2(x, y, z); //Creates second vector
			cout << "Vector: " << V2 << "initialized!" << endl;
            out << "Vector: " << V2 << "initialized!" << endl;
            if (V1.perpendicular_vectors(V2)) //Prints result based on entered vector points
            {
                cout << "\tVectors " << V1 << " and " << V2 << " are perpendicular" << endl << endl;
                out << "\tVectors " << V1 << " and " << V2 << " are perpendicular" << endl << endl;

            }
            else
            {
                cout << "\tVectors " << V1 << " and " << V2 << " are NOT perpendicular" << endl << endl;
                out << "\tVectors " << V1 << " and " << V2 << " are NOT perpendicular" << endl << endl;
            }
            break;
		}
		case 6: //Sums two vectors
		{
			input >> x;
			input >> y;
			input >> z;
			Vector V2(x, y, z); //Creates second vector
			cout << "Vector: " << V2 << "initialized!" << endl;
            out << "Vector: " << V2 << "initialized!" << endl;

            Vector V3 = V1 + V2; //Sums the vectors
			cout << "Vector " << V1 << " + Vector " << V2 << " is equal to Vector " << V3 << endl << endl; //Prints result
			out << "Vector " << V1 << " + Vector " << V2 << " is equal to Vector " << V3 << endl << endl; //Prints result
			break;
		}
		case 7: //Multiplies vector with a number
		{
			float x;
			input >> x; //Reads the number
			Vector V2 = V1 * x; //Multiplies them
			cout << "Vector " << V1 << "multiplied with " << x << " is equal to Vector " << V2 << endl << endl; //Prints the result
			out << "Vector " << V1 << "multiplied with " << x << " is equal to Vector " << V2 << endl << endl; //Prints the result
			break;
		}
		case 8: //Scalar multiplication of two vectors
		{
			input >> x;
			input >> y;
			input >> z;
			Vector V2(x, y, z); //Creates second vector
			cout << "Vector: " << V2 << "initialized!" << endl;
            out << "Vector: " << V2 << "initialized!" << endl;

            float result = V1 * V2; //Multiplies the vectors
			cout << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to " << result << endl << endl; //Prints the result
			out << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to " << result << endl << endl; //Prints the result
			break;
		}
		case 9:
		{
			input >> x;
			input >> y;
			input >> z;
			Vector V2(x, y, z); //Creates second vector
			cout << "Vector: " << V2 << "initialized!" << endl;
            out << "Vector: " << V2 << "initialized!" << endl;

            Vector V3 = V1 ^ V2; //Multiplies the vectors scalarly
			cout << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to Vector " << V3 << endl << endl; //Prints the result
			out << "Vector " << V1 << " scalar multiplied to Vector " << V2 << " is equal to Vector " << V3 << endl << endl; //Prints the result
			break;
		}
		case 10: //Mixed multiplication of vectors
		{
			input >> x;
			input >> y;
			input >> z;
			Vector V2(x, y, z); //Creates second vector
			cout << "Vector: " << V2 << "initialized!" << endl;
            out << "Vector: " << V2 << "initialized!" << endl;

			input >> x;
			input >> y;
			input >> z;
			Vector V3(x, y, z); //Creates third vector
            cout << "Vector: " << V3 << "initialized!" << endl;
            out << "Vector: " << V3 << "initialized!" << endl;

            float result = V1.operator()(V2, V3);
			cout << "Mixed multiplication of Vector " << V1 << " Vector " << V2 << " and Vector " << V3 << " is equal to " << result << endl << endl; //Does the operation and prints the result
			out << "Mixed multiplication of Vector " << V1 << " Vector " << V2 << " and Vector " << V3 << " is equal to " << result << endl << endl; //Does the operation and prints the result
			break;
		}
		default:
        {
            cout << "Invalid choice!" << endl;
            out << "Invalid choice!" << endl;
			break;
        } //If input is invalid
    }
		other_operation = false; //Breaks the while loop
		cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
		char c;
		cin >> c;
		if (c == 'Y') {
			other_object = true; //Returns to main menu
		}
		else {
			other_object = false; //Exits program
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}

}


void Menu::lineMenu() //Line menu
{
	if (readData(choose)) //If user wants to read from console
	{
		float x, y, z;
		int choice;

		out.open("output.txt");
		cout << "\tEnter values for x, y, z to define Point P1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P1(x, y, z); //Creating creating point P1
		cout << "Point P1:" << P1 << "initialized!" << endl;
		out << "Point P1:" << P1 << "initialized!" << endl;

		cout << "\tEnter values for x, y, z to define Vector V1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Vector V1(x, y, z); //Creating vector V1
		cout << "Vector :" << V1 << "initialized!" << endl;
		out << "Vector :" << V1 << "initialized!" << endl;

		Line l1(P1, V1);
		cout << "Line l1 is defined by Point " << P1 << " and Vector " << V1 << endl << endl;
		out << "Line l1 is defined by Point " << P1 << " and Vector " << V1 << endl << endl;

		while (other_operation)
		{

			cout << "************************** L I N E **************************" << endl;
			cout << endl << endl << endl;
			cout << "PLEASE CHOOSE OPERATION FROM BELOW OR PRESS 0 TO GO BACK TO MAIN MENU" << endl << endl;
			cout << "\t1.Find line direction" << endl << endl;
			cout << "\t2.Find normal vector" << endl << endl;
			cout << "\t3.Find angle between two lines" << endl << endl;
			cout << "\t4.Check if a Point lays on the Line" << endl << endl;
			cout << "\t5.Check if a Line is parallel to other Line" << endl << endl;
			cout << "\t6.Check if a Line coincide to other Line" << endl << endl;
			cout << "\t7.Check if a Line crosses other Line" << endl << endl;
			cout << "\t8.Check if a Line is crossed with other Line" << endl << endl;
			cout << "\t9.Check if a Line is perpendicular with other Line" << endl << endl;

			cin >> choice; //Reading user's choice

			if (choice == 0) //If he wants to go back to main menu - break the while loop
			{
				break;
			}

			switch (choice) { //Checks what the user's input is and does the following action
			case 1: //Finds line direction
			{
				cout << "Line direction is Vector " << l1.find_direction() << endl; //Prints line direction
				out << "Line direction is Vector " << l1.find_direction() << endl; //Prints line direction
				break;
			}
			case 2: //Finds normal vector
			{
				cout << "Normal vector is " << l1.normal_vector() << endl; //Prints normal vector
				out << "Normal vector is " << l1.normal_vector() << endl; //Prints normal vector
				break;
			}
			case 3: //Finds angle between two lines
			{
				cout << "\tEnter values for x, y, z to define Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P2(x, y, z); //Creates point 2
				cout << "Point: " << P1 << "initialized!" << endl;
				out << "Point: " << P1 << "initialized!" << endl;

				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates vector 2
                cout << "Vector: " << V2 << "initialized!" << endl;
                out << "Vector: " << V2 << "initialized!" << endl;

                Line l2(P2, V2); //Creates line 2
				cout << "Angle between Line " << P1 << " " << V1 << " and Line "
					<< P1 << "" << V1 << " is " << l1.Angel_between_two_lines(l2)
					<< " radians" << endl; //Prints angle between the lines
				out << "Angle between Line " << P1 << " " << V1 << " and Line "
					<< P1 << "" << V1 << " is " << l1.Angel_between_two_lines(l2)
					<< " radians" << endl;
                break;
			}
			case 4: //Checks if a point lays on the line
			{
				cout << "\tEnter values for x, y, z to define Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P2(x, y, z); //Creates point 2
				cout << "Point: " << P1 << "initialized!" << endl;
				out << "Point: " << P1 << "initialized!" << endl;

                if (l1.operator+(P2)) //Prints result
                {
                    cout << "Point " << P2 << " LAYS on Line" << endl;
                    out << "Point " << P2 << " LAYS on Line" << endl;
                }
				else
                {
                    cout << "Point " << P2 << " DOES NOT LAY on Line" << endl;
					out << "Point " << P2 << " DOES NOT LAY on Line" << endl;
                }
				break;
			}
			case 5: //Checks if two lines are parallel
			{
				cout << "\tEnter values for x, y, z to define Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2 << "initialized!" << endl;


				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

                Line l2(P2, V2); //Creates second line

				if (l1.operator||(l2)) //Prints result
                {
                    cout << "Line " << P1 << ", " << V1 << " is prallel to Line " << P2 << ", " << V2 << endl;
                    out << "Line " << P1 << ", " << V1 << " is prallel to Line " << P2 << ", " << V2 << endl;
                }
				else
                {
                    cout << "Line " << P1 << ", " << V1 << " is not parallel to Line " << P2 << ", " << V2 << endl;
                    out << "Line " << P1 << ", " << V1 << " is not parallel to Line " << P2 << ", " << V2 << endl;
                }
				break;
			}
			case 6: //Checks if two lines coincide
			{
				cout << "\tEnter values for x, y, z to define Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2<< "initialized!" << endl;

				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

				Line l2(P2, V2); //Creates second line

				if (l1.operator==(l2))//Prints result based on outcome
				{
                    cout << "Line " << P1 << ", " << V1 << " coincides with Line " << P2 << ", " << V2 << endl;
					out << "Line " << P1 << ", " << V1 << " coincides with Line " << P2 << ", " << V2 << endl;
				}
				else
                {
					cout << "Line " << P1 << ", " << V1 << " does not coincide with Line " << P2 << ", " << V2 << endl;
					out << "Line " << P1 << ", " << V1 << " does not coincide with Line " << P2 << ", " << V2 << endl;
                }
				break;
			}
			case 7: //Checks if two lines cross eachother
			{
				cout << "\tEnter values for x, y, z to define Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2<< "initialized!" << endl;

				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

				Line l2(P2, V2); //Creates second line

				if (l1.operator||(l2)) //Prints result based on outcome
				{
                    cout << "Line " << P1 << " " << V1 << " does not cross Line " << P2 << " " << V2 << endl;
                    out << "Line " << P1 << " " << V1 << " does not cross Line " << P2 << " " << V2 << endl;
				}
				else
                {
                    cout << "Line " << P1 << " " << V1 << " crosses Line " << P2 << " " << V2 << endl;
                    out << "Line " << P1 << " " << V1 << " crosses Line " << P2 << " " << V2 << endl;
                }
				break;
			}
			case 9: // Checks if two lines are perpendicular
			{
				cout << "\tEnter values for x, y, z to define Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2<< "initialized!" << endl;

				cout << "\tEnter values for x, y, z to define Vector V2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

				Line l2(P2, V2); //Creates second line

				if (V1.perpendicular_vectors(V2)) //Prints result based on outcome
				{
                    cout << "Line " << l1 << " IS perpendicular to Line " << l2 << endl;
  					out << "Line " << l1 << " IS perpendicular to Line " << l2 << endl;
				}
				else
                {
					cout << "Line " << l1 << " IS NOT perpendicular to Line " << l2 << endl;
					out << "Line " << l1 << " IS NOT perpendicular to Line " << l2 << endl;
                }
				break;
			}
			default: //If input < 1 || input > 9
            {
			   cout << "Invalid input!" << endl;
			   out << "Invalid input!" << endl;
			   break;
            }
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //Goes back to line menu
			else {
				other_operation = false; //Breaks the while loop
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y')
					other_object = true; //Goes back to main menu
				else
					other_object = false; //Exits program
			}
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}
	else // If user wants to read from a file
	{
		ifstream input("Line.txt"); //Opens Line.txt file
		if (!input) //If input is invalid
		{
			cerr << "File cannot be opened!" << endl;
		}

		out.open("output.txt");

		input.seekg(1, ios::beg);
		float x, y, z;
		input >> x;
		input >> y;
		input >> z;
		Point P1(x, y, z); //Creates point

		input >> x;
		input >> y;
		input >> z;
		Vector V1(x, y, z); //Creates vector

		Line l1(P1, V1); //Creates line one
		cout << "Line L1 initialized: " << P1 << "" << V1 << endl << endl; //Prints line one

		int choice;
		input >> choice;

		switch (choice) { //Based on the file input does the needed operation
		case 1: //Finds line direction
			{
				cout << "Line direction is Vector " << l1.find_direction() << endl; //Prints line direction
				out << "Line direction is Vector " << l1.find_direction() << endl; //Prints line direction
				break;
			}
			case 2: //Finds normal vector
			{
				cout << "Normal vector is " << l1.normal_vector() << endl; //Prints normal vector
				out << "Normal vector is " << l1.normal_vector() << endl; //Prints normal vector
				break;
			}
			case 3: //Finds angle between two lines
			{
				input >> x;
				input >> y;
				input >> z;
				Point P2(x, y, z); //Creates point 2
				cout << "Point: " << P1 << "initialized!" << endl;
				out << "Point: " << P1 << "initialized!" << endl;

				input >> x;
				input >> y;
				input >> z;
				Vector V2(x, y, z); //Creates vector 2
                cout << "Vector: " << V2 << "initialized!" << endl;
                out << "Vector: " << V2 << "initialized!" << endl;

                Line l2(P2, V2); //Creates line 2
				cout << "Angle between Line " << P1 << " " << V1 << " and Line "
					<< P1 << "" << V1 << " is " << l1.Angel_between_two_lines(l2)
					<< " radians" << endl; //Prints angle between the lines
				out << "Angle between Line " << P1 << " " << V1 << " and Line "
					<< P1 << "" << V1 << " is " << l1.Angel_between_two_lines(l2)
					<< " radians" << endl;
                break;
			}
			case 4: //Checks if a point lays on the line
			{
				input >> x;
				input >> y;
				input >> z;
				Point P2(x, y, z); //Creates point 2
				cout << "Point: " << P1 << "initialized!" << endl;
				out << "Point: " << P1 << "initialized!" << endl;

                if (l1.operator+(P2)) //Prints result
                {
                    cout << "Point " << P2 << " LAYS on Line" << endl;
                    out << "Point " << P2 << " LAYS on Line" << endl;
                }
				else
                {
                    cout << "Point " << P2 << " DOES NOT LAY on Line" << endl;
					out << "Point " << P2 << " DOES NOT LAY on Line" << endl;
                }
				break;
			}
			case 5: //Checks if two lines are parallel
			{
				input >> x;
				input >> y;
				input >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2 << "initialized!" << endl;


				input >> x;
				input >> y;
				input >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

                Line l2(P2, V2); //Creates second line

				if (l1.operator||(l2)) //Prints result
                {
                    cout << "Line " << P1 << ", " << V1 << " is prallel to Line " << P2 << ", " << V2 << endl;
                    out << "Line " << P1 << ", " << V1 << " is prallel to Line " << P2 << ", " << V2 << endl;
                }
				else
                {
                    cout << "Line " << P1 << ", " << V1 << " is not parallel to Line " << P2 << ", " << V2 << endl;
                    out << "Line " << P1 << ", " << V1 << " is not parallel to Line " << P2 << ", " << V2 << endl;
                }
				break;
			}
			case 6: //Checks if two lines coincide
			{
				input >> x;
				input >> y;
				input >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2<< "initialized!" << endl;

				input >> x;
				input >> y;
				input >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

				Line l2(P2, V2); //Creates second line

				if (l1.operator==(l2))//Prints result based on outcome
				{
                    cout << "Line " << P1 << ", " << V1 << " coincides with Line " << P2 << ", " << V2 << endl;
					out << "Line " << P1 << ", " << V1 << " coincides with Line " << P2 << ", " << V2 << endl;
				}
				else
                {
					cout << "Line " << P1 << ", " << V1 << " does not coincide with Line " << P2 << ", " << V2 << endl;
					out << "Line " << P1 << ", " << V1 << " does not coincide with Line " << P2 << ", " << V2 << endl;
                }
				break;
			}
			case 7: //Checks if two lines cross eachother
			{
				input >> x;
				input >> y;
				input >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2<< "initialized!" << endl;

				input >> x;
				input >> y;
				input >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

				Line l2(P2, V2); //Creates second line

				if (l1.operator||(l2)) //Prints result based on outcome
				{
                    cout << "Line " << P1 << " " << V1 << " does not cross Line " << P2 << " " << V2 << endl;
                    out << "Line " << P1 << " " << V1 << " does not cross Line " << P2 << " " << V2 << endl;
				}
				else
                {
                    cout << "Line " << P1 << " " << V1 << " crosses Line " << P2 << " " << V2 << endl;
                    out << "Line " << P1 << " " << V1 << " crosses Line " << P2 << " " << V2 << endl;
                }
				break;
			}
			case 9: // Checks if two lines are perpendicular
			{
				input >> x;
				input >> y;
				input >> z;
				Point P2(x, y, z); //Creates point
				cout << "Point: " << P2 << "initialized!" << endl;
				out << "Point: " << P2<< "initialized!" << endl;

				input >> x;
				input >> y;
				input >> z;
				Vector V2(x, y, z); //Creates vector
				cout << "Vector: " << V2 << "initialized!" << endl;
				out << "Vector: " << V2 << "initialized!" << endl;

				Line l2(P2, V2); //Creates second line

				if (V1.perpendicular_vectors(V2)) //Prints result based on outcome
				{
                    cout << "Line " << l1 << " IS perpendicular to Line " << l2 << endl;
  					out << "Line " << l1 << " IS perpendicular to Line " << l2 << endl;
				}
				else
                {
					cout << "Line " << l1 << " IS NOT perpendicular to Line " << l2 << endl;
					out << "Line " << l1 << " IS NOT perpendicular to Line " << l2 << endl;
                }
				break;
			}
			default: //If input < 1 || input > 9
            {
			   cout << "Invalid input!" << endl;
			   out << "Invalid input!" << endl;
			   break;
            }
        }
		other_operation = false; //Breaks the while loop
		cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
		char c;
		cin >> c;
		if (c == 'Y') {
			other_object = true; //Goes to main menu
		}
		else {
			other_object = false; //Exits program
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}
}

void Menu::segmentMenu()
{
	if (readData(choose)) //If user wants to enter from the console
	{
	    out.open("output.txt");
		float x, y, z;
		int choice;
		cout << "\tEnter values for x, y, z to define the starting Point P1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P1(x, y, z); //Creates first point
		cout << "Point: " << P1 << "initialized!" << endl;
		out << "Point: " << P1 << "initialized!" << endl;

		cout << "\tEnter values for x, y, z to define the ending Point P2" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P2(x, y, z); //Creates second point
		cout << "Point: " << P2 << "initialized!" << endl;
		out << "Point: " << P2 << "initialized!" << endl;

		Segment S(P1, P2); //Creates the segment

		bool flag = true;

		while (other_operation)
		{
			cout << "**************************  S E G M E N T  **************************" << endl;
			cout << endl << endl << endl;
			cout << "PLEASE, CHOOSE OPERATION FROM BELOW OR PRESS 0 TO GO BACK TO MAIN MENU" << endl << endl;
			cout << "\t1.Find the length of the Segment " << endl << endl;
			cout << "\t2.Find the middle of the Segment" << endl << endl;
			cout << "\t3.Does a Point lays on the Segment" << endl << endl;

			cin >> choice;

			if (choice == 0) //If input is 0, goes back to main menu
			{
				break;
			}

			switch (choice) { //Doing the following operation based on the users input
			case 1: //Finds the length of the segment
			{
				cout << "The length of the Segment is " << S.find_length_of_Segment() << endl; //Prints result
				out << "The length of the Segment is " << S.find_length_of_Segment() << endl; //Prints result
				break;
			}
			case 2: //Finds the middle of the segment
			{
				cout << "The middle of the Segment is " << S.find_Mid_point_of_Segment() << endl; //Pritns result
				out << "The middle of the Segment is " << S.find_Mid_point_of_Segment() << endl; //Pritns result
				break;
			}
			case 3: //Checks if a point lays on the segment
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P3(x, y, z); //Creates the point
				cout << "Point: " << P3 << "initialized!" << endl;
				out << "Point: " << P3 << "initialized!" << endl;

 				if (S.operator==(P3)) //Checks if the point lays on the segment and prints the result
				{
				    cout << "Point " << P3 << " LAYS on the Segment " << endl;
				    out << "Point " << P3 << " LAYS on the Segment " << endl;
				}
				else
                {
                    cout << "Point " << P3 << " DOES NOT LAY on the Segment " << endl;
					out << "Point " << P3 << " DOES NOT LAY on the Segment " << endl;
                }
				break;
			}
			default:
            {
			   cout << "Invalid choice!" << endl;
			   out << "Invalid choice!" << endl;
			   break;
            } //If input is invalid
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //Goes back to segment menu
			else {
				other_operation = false;
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y')
					other_object = true; //Goes back to main menu
				else
					other_object = false; //Exits program
			}
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}
	else //User wants to read from a file
	{
		ifstream input("Segment.txt"); //Opens segment.txt and reads the information from it
		if (!input) //If file doesnt exist
		{
			cerr << "File cannot be opened!" << endl;
		}
		input.seekg(1, ios::beg);

		out.open("output.txt");

		float x, y, z;
		input >> x;
		input >> y;
		input >> z;
		Point P1(x, y, z); //Creates first point
		cout << "Point: " << P1 << "initialized!" << endl;
		out << "Point: " << P1 << "initialized!" << endl;

		input >> x;
		input >> y;
		input >> z;
		Point P2(x, y, z); //Creates second point
		cout << "Point: " << P2 << "initialized!" << endl;
		out << "Point: " << P2 << "initialized!" << endl;

		Segment S(P1, P2); //Creates the segment
		cout << "Segment S1 initialized: " << P1 << "" << P2 << endl << endl;

		int choice;
		input >> choice;

		switch (choice) { //Reads from the file which operation should be done and does it if it is valid
		case 1: //Finds the length of the segment
		{
			cout << "The length of the Segment is " << S.find_length_of_Segment() << endl; //Prints the result
			out << "The length of the Segment is " << S.find_length_of_Segment() << endl; //Prints the result
			break;
		}
		case 2: //Finds the middle of the segment
		{
			cout << "The middle of the Segment is " << S.find_Mid_point_of_Segment() << endl; //Prints the result
			out << "The middle of the Segment is " << S.find_Mid_point_of_Segment() << endl; //Prints the result
			break;
		}
		case 3: // Checks if a point lays on the segment
		{
			input >> x;
			input >> y;
			input >> z;
			Point P3(x, y, z); //Creates the point
			cout << "Point: " << P3 << "initialized!" << endl;
			out << "Point: " << P3 << "initialized!" << endl;

			if (S.operator==(P3)) //Checks if the point lays on the segment and prints the result
			{
			    cout << "Point " << P3 << " LAYS on the Segment " << endl;
				out << "Point " << P3 << " LAYS on the Segment " << endl;
			}
			else
            {
                cout << "Point " << P3 << " DOES NOT LAY on the Segment " << endl;
				out << "Point " << P3 << " DOES NOT LAY on the Segment " << endl;
            }
			break;
		}
		default: //If input is invalid
        {
		   cout << "Invalid choice!" << endl;
		   out << "Invalid choice!" << endl;
		   break;
        }
    }
		other_operation = false; //Breaks the while loop
		cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
		char c;
		cin >> c;
		if (c == 'Y') {
			other_object = true; //Goes to main menu
		}
		else {
			other_object = false; //Exits the program
		}
	}
	out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
}


void Menu::triangleMenu() //Triangle menu
{
	if (readData(choose)) //If user wants to read from the console
	{
	    out.open("output.txt");
		float x, y, z;
		int choice;
		cout << "\tEnter values for x, y, z to define the starting Point P1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P1(x, y, z); //Creates first point
		cout << "Point: " << P1 << "initialized!" << endl;
        out << "Point: " << P1 << "initialized!" << endl;

		cout << "\tEnter values for x, y, z to define the ending Point P2" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P2(x, y, z); //Creates second point
		cout << "Point: " << P2 << "initialized!" << endl;
        out << "Point: " << P2 << "initialized!" << endl;

        cout << "\tEnter values for x, y, z to define the ending Point P3" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P3(x, y, z); //Creates third point
		cout << "Point: " << P3 << "initialized!" << endl;
        out << "Point: " << P3 << "initialized!" << endl;

		Triangle T(P1, P2, P3); //Creates the triangle

		while (other_operation)
		{
			cout << "**************************  T R I A N G L E  **************************" << endl;
			cout << endl << endl << endl;
			cout << "PLEASE, CHOOSE OPERATION FROM BELOW OR PRESS 0 TO GO BACK TO MAIN MENU" << endl << endl;
			cout << "\t1.Define the type of the Triangle " << endl << endl;
			cout << "\t2.Find the Area of the Triangle" << endl << endl;
			cout << "\t3.Find the Perimeter of the Triangle" << endl << endl;
			cout << "\t4.Find the Medicenter of the Triangle" << endl << endl;
			cout << "\t5.Check if Point lays within Triangle and is inside Triangle" << endl << endl;
			cout << "\t6.Check if Point lays within Triangle and is outside Triangle" << endl << endl;
			cout << "\t7.Check if Point lays on a Segment of the Triangle" << endl << endl;

			cin >> choice;

			if (choice == 0) //If input is 0, goes back to the main menu
			{
				break;
			}

			switch (choice) { //Checks the input and does the following operation
			case 1: //Defines the type of the triangle
			{
				T.define_Type(); //Prints result
				break;
			}
			case 2: //Finds the area of the triangle
			{
				cout << "Area of the Triangle is " << T.find_Area() << endl; //Prints result
				out << "Area of the Triangle is " << T.find_Area() << endl; //Prints result
				break;
			}
			case 3: //Finds the perimeter of the triangle
			{
				cout << "Perimeter of the Triangle is " << T.find_Perimeter() << endl; //Prints result
				out << "Perimeter of the Triangle is " << T.find_Perimeter() << endl; //Prints result
				break;
			}
			case 4: //Find the medicenter of the triangle
			{
				cout << "Medicenter of the Triangle is " << T.find_Medicenter() << endl; //Prints the result
				out << "Medicenter of the Triangle is " << T.find_Medicenter() << endl; //Prints the result
				break;
			}
			case 5: // Checks if a point lays on the surface of the triangle and is posessed inside it
			{
				cout << "\tEnter values for x, y, z to define the Point P" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (T.operator<(P)) //Checks if the point lays on the surface of the triangle and is possessed inside it and then prints the result
				{
                    cout << "Point " << P << " LAYS on the surface of the Triangle and is possessed INSIDE of it" << endl;
 					out << "Point " << P << " LAYS on the surface of the Triangle and is possessed INSIDE of it" << endl;
				}
				else{
					cout << "Point " << P << " DOES NOT LAY on the surface of the Triangle or" << endl;
					cout << " is not possessed INSIDE of it or both" << endl;
					out << "Point " << P << " DOES NOT LAY on the surface of the Triangle or" << endl;
					out << " is not possessed INSIDE of it or both" << endl;
				}
				break;
			}
			case 6: // Checks if a point lays on the surface of the triangle and is possessed outside of it
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (T.operator>(P)) //Checks if the point lays on the surface and is possessed outside of it and then prints the result
                {
                    cout << "Point " << P << " LAYS on the surface of the Triangle and is possessed OUTSIDE of it";
                    out << "Point " << P << " LAYS on the surface of the Triangle and is possessed OUTSIDE of it";
                }
				else{
					cout << "Point " << P << "DOES NOT LAY on the surface of the Triangle or" << endl;
					cout << " is not possessed OUTSIDE of it or both" << endl;
					out << "Point " << P << "DOES NOT LAY on the surface of the Triangle or" << endl;
					out << " is not possessed OUTSIDE of it or both" << endl;
				}
				break;
			}
			case 7: // Checks if a point lays on a segment of the triangle
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (T.operator==(P)) //Checks if the point lays on the segment of the triangle and prints the result
                {
                    cout << "Point " << P << " lays on a Segment of the Triangle" << endl;
                    out << "Point " << P << " lays on a Segment of the Triangle" << endl;
                }
				else
                {
                    cout << "Point " << P << " doesnt lay on a Segment of the Triangle" << endl;
					out << "Point " << P << " doesnt lay on a Segment of the Triangle" << endl;
                }
				break;
			}
			default: //If input is invalid
            {
                cout << "Invalid choice!" << endl;
                out << "Invalid choice!" << endl;
				break;
            }
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //Goes to triangle menu
			else {
				other_operation = false; //Breaks the while loop
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y')
					other_object = true; //Goes back to the main menu
				else
					other_object = false; //Exits program
			}
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
    }
	else //User wants to read from a file
	{
		ifstream input("Triangle.txt"); //Opens file Triangle.txt
		if (!input) //If file is invalid
		{
			cerr << "File cannot be opened!" << endl;
		}
		input.seekg(1, ios::beg);

		out.open("output.txt");

		float x, y, z;
		input >> x;
		input >> y;
		input >> z;
		Point P1(x, y, z); //Creates first point
		cout << "Point: " << P1 << "initialized!" << endl;
        out << "Point: " << P1 << "initialized!" << endl;

		input >> x;
		input >> y;
		input >> z;
		Point P2(x, y, z); //Creates second point
		cout << "Point: " << P2 << "initialized!" << endl;
        out << "Point: " << P2 << "initialized!" << endl;

		input >> x;
		input >> y;
		input >> z;
		Point P3(x, y, z); //Creates third point
		cout << "Point: " << P3 << "initialized!" << endl;
        out << "Point: " << P3 << "initialized!" << endl;

		Triangle T(P1, P2, P3); //Creates the triangle
		cout << "Triangle T initialized with: " << P1 << " " << P2 << "" << P3 << endl << endl;
		out << "Triangle T initialized with: " << P1 << " " << P2 << "" << P3 << endl << endl;

		int choice;
		input >> choice;

		switch (choice) { //Checks the file input and does the following action
		case 1: //Defines the type of the triangle
			{
				T.define_Type(); //Prints result
				break;
			}
			case 2: //Finds the area of the triangle
			{
				cout << "Area of the Triangle is " << T.find_Area() << endl; //Prints result
				out << "Area of the Triangle is " << T.find_Area() << endl; //Prints result
				break;
			}
			case 3: //Finds the perimeter of the triangle
			{
				cout << "Perimeter of the Triangle is " << T.find_Perimeter() << endl; //Prints result
				out << "Perimeter of the Triangle is " << T.find_Perimeter() << endl; //Prints result
				break;
			}
			case 4: //Find the medicenter of the triangle
			{
				cout << "Medicenter of the Triangle is " << T.find_Medicenter() << endl; //Prints the result
				out << "Medicenter of the Triangle is " << T.find_Medicenter() << endl; //Prints the result
				break;
			}
			case 5: // Checks if a point lays on the surface of the triangle and is posessed inside it
			{
				cout << "\tEnter values for x, y, z to define the Point P" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (T.operator<(P)) //Checks if the point lays on the surface of the triangle and is possessed inside it and then prints the result
				{
                    cout << "Point " << P << " LAYS on the surface of the Triangle and is possessed INSIDE of it" << endl;
 					out << "Point " << P << " LAYS on the surface of the Triangle and is possessed INSIDE of it" << endl;
				}
				else{
					cout << "Point " << P << " DOES NOT LAY on the surface of the Triangle or" << endl;
					cout << " is not possessed INSIDE of it or both" << endl;
					out << "Point " << P << " DOES NOT LAY on the surface of the Triangle or" << endl;
					out << " is not possessed INSIDE of it or both" << endl;
				}
				break;
			}
			case 6: // Checks if a point lays on the surface of the triangle and is possessed outside of it
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (T.operator>(P)) //Checks if the point lays on the surface and is possessed outside of it and then prints the result
                {
                    cout << "Point " << P << " LAYS on the surface of the Triangle and is possessed OUTSIDE of it";
                    out << "Point " << P << " LAYS on the surface of the Triangle and is possessed OUTSIDE of it";
                }
				else{
					cout << "Point " << P << "DOES NOT LAY on the surface of the Triangle or" << endl;
					cout << " is not possessed OUTSIDE of it or both" << endl;
					out << "Point " << P << "DOES NOT LAY on the surface of the Triangle or" << endl;
					out << " is not possessed OUTSIDE of it or both" << endl;
				}
				break;
			}
			case 7: // Checks if a point lays on a segment of the triangle
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (T.operator==(P)) //Checks if the point lays on the segment of the triangle and prints the result
                {
                    cout << "Point " << P << " lays on a Segment of the Triangle" << endl;
                    out << "Point " << P << " lays on a Segment of the Triangle" << endl;
                }
				else
                {
                    cout << "Point " << P << " doesnt lay on a Segment of the Triangle" << endl;
					out << "Point " << P << " doesnt lay on a Segment of the Triangle" << endl;
                }
				break;
			}
			default: //If input is invalid
            {
                cout << "Invalid choice!" << endl;
                out << "Invalid choice!" << endl;
				break;
            }
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //Goes to triangle menu
			else {
				other_operation = false; //Breaks the while loop
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y')
					other_object = true; //Goes back to the main menu
				else
					other_object = false; //Exits program
			}
			out.close();
			if(other_object == false)
            {
                cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
                char write;
                cin >> write;
                if(write != '1')
                {
                    remove("output.txt");
                }
            }
		}
}

void Menu::tetrahedronMenu() //Tetrahedron menu
{
	if (readData(choose)) //User wants to read from the console
	{
	    out.open("output.txt");
		float x, y, z;
		int choice;
		cout << "\tEnter values for x, y, z to define the starting Point P1" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P1(x, y, z); //Creates the first point
		cout << "Point: " << P1 << "initialized!" << endl;
        out << "Point: " << P1 << "initialized!" << endl;

		cout << "\tEnter values for x, y, z to define the ending Point P2" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P2(x, y, z); //Creates the second point
		cout << "Point: " << P2 << "initialized!" << endl;
        out << "Point: " << P2 << "initialized!" << endl;

		cout << "\tEnter values for x, y, z to define the ending Point P3" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P3(x, y, z); //Creates the third point
		cout << "Point: " << P3 << "initialized!" << endl;
        out << "Point: " << P3 << "initialized!" << endl;

		cout << "\tEnter values for x, y, z to define the ending Point P4" << endl;
		cout << "\tx = ";
		cin >> x;
		cout << "\ty = ";
		cin >> y;
		cout << "\tz = ";
		cin >> z;
		Point P4(x, y, z); //Creates the fourth point
		cout << "Point: " << P4 << "initialized!" << endl;
        out << "Point: " << P4 << "initialized!" << endl;

		Tetrahedron Te(P1, P2, P3, P4); //Creates the tetrahedron

		while (other_operation)
		{
			cout << "**************************  T E T R A H E D R O N  **************************" << endl;
			cout << endl << endl << endl;
			cout << "PLEASE, CHOOSE AN OPERATION FROM BELOW OR PRESS 0 TO GO BACK TO MAIN MENU" << endl << endl;
			cout << "\t1.Check if Tetrahedron is regular" << endl << endl;
			cout << "\t2.Check if Tetrahedron is orthogonal" << endl << endl;
			cout << "\t3.Find surface of Tetrahedron" << endl << endl;
			cout << "\t4.Find volume of Tetrahedron" << endl << endl;
			cout << "\t5.Check if Point lays within Tetrahedron and is inside Tetrahedron" << endl << endl;
			cout << "\t6.Check if Point lays within Tetrahedron and is outside Tetrahedron" << endl << endl;
			cout << "\t7.Check if Point lays on a Segment of the Tetrahedron" << endl << endl;


			cin >> choice;

			if (choice == 0) //Chcks if user wants to go back to main menu
			{
				break;
			}

			switch (choice) { //Takes the user's input and does the following operation
			case 1: //Checks if the tetrahedron is regular
			{
				if (Te.check_if_Regular())
                {
                    cout << "Tetrahedron is regular" << endl; //Prints the result
   					out << "Tetrahedron is regular" << endl; //Prints the result
                }
				else
                {
 					cout << "Tetrahedron is irregular" << endl; //Prints the result
					out << "Tetrahedron is irregular" << endl; //Prints the result
                }
				break;
			}
			case 2: //Checks if the tetrahedron is orthogonal
			{
				if (Te.check_if_Ortogonal())
                {
 					cout << "Tetrahedron is orthogonal" << endl; //Prints the result
					out << "Tetrahedron is orthogonal" << endl; //Prints the result
                }
				else
                {
 					cout << "Tetrahedron is not orthogonal" << endl; //Prints the result
					out << "Tetrahedron is not orthogonal" << endl; //Prints the result
                }
				break;
			}
			case 3: //Finds the surface of the tetrahedron
			{
				cout << "Surface of the Tetrahedron is " << Te.get_Surface() << endl; //Prints the result
				out << "Surface of the Tetrahedron is " << Te.get_Surface() << endl; //Prints the result
				break;
			}
			case 4: //Finds the volume of the tetrahedron
			{
				cout << "Volume of the Tetrahedron is " << Te.get_Volume() << endl; //Prints the result
				out << "Volume of the Tetrahedron is " << Te.get_Volume() << endl; //Prints the result
				break;
			}
			case 5: //Checks if a point lays on the surface of the tetrahedron and is possessed inside of it
			{
				cout << "\tEnter values for x, y, z to define the Point P" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (Te.operator<(P))
                {
                    cout << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed INSIDE it" << endl;  //Prints the result
                    out << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed INSIDE it" << endl;  //Prints the result
                }
				else{
					cout << "Point " << P << " DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					cout << "  is not possessed INSIDE it or both" << endl;  //Prints the result
					out << "Point " << P << " DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					out << "  is not possessed INSIDE it or both" << endl;
				}
				break;
			}
			case 6: //Checks if a point lays on the surface of the tetrahedron and is possessed outside of it
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (Te.operator>(P))
                {
                    cout << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed OUTSIDE of it" << endl;  //Prints the result
					out << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed OUTSIDE of it" << endl;  //Prints the result
                }
				else{
					cout << "Point " << P << "DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					cout << " is not possessed OUTSIDE of it or both" << endl;  //Prints the result
					out << "Point " << P << "DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					out << " is not possessed OUTSIDE of it or both" << endl;
				}
				break;
			}
			case 7: //Checks if a point lays on the segment of the tetrahedron
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (Te.operator==(P))
                {
                    cout << "Point " << P << " lays on a Segment of the Tetrahedron" << endl; //Prints the result
                    out << "Point " << P << " lays on a Segment of the Tetrahedron" << endl; //Prints the result
                }
				else
                {
 					cout << "Point " << P << " does not lay on a Segment of the Triangle" << endl; //Prints the result
                    out << "Point " << P << " does not lay on a Segment of the Triangle" << endl; //Prints the result
                }
				break;
			}
			default: //If input is invalid
            {
                cout << "Invalid choice!" << endl;
                out << "Invalid choice!" << endl;
                break;
            }
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //Goes back to tetrahedron menu
			else{
				other_operation = false; //Breaks the while loop
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y')
					other_object = true; //Goes back to the main menu
				else
					other_object = false; //Exit program
			}
		}
		out.close();
		if(other_object == false)
        {
            cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
            char write;
            cin >> write;
            if(write != '1')
            {
                remove("output.txt");
            }
        }
	}
	else //User wants to read from a file
	{
		ifstream input("Tetrahedron.txt"); //Opens Tetrahedron.txt
		if (!input) //Checks if file is valid
		{
			cerr << "File cannot be opened!" << endl;
		}
		input.seekg(1, ios::beg);

		out.open("output.txt");
		float x, y, z;
		input >> x;
		input >> y;
		input >> z;
		Point P1(x, y, z); //Creates first point
		cout << "Point: " << P1 << "initialized!" << endl;
        out << "Point: " << P1 << "initialized!" << endl;

		input >> x;
		input >> y;
		input >> z;
		Point P2(x, y, z); //Creates second point
		cout << "Point: " << P2 << "initialized!" << endl;
        out << "Point: " << P2 << "initialized!" << endl;

		input >> x;
		input >> y;
		input >> z;
		Point P3(x, y, z); //Creates third point
		cout << "Point: " << P3 << "initialized!" << endl;
        out << "Point: " << P3 << "initialized!" << endl;

		input >> x;
		input >> y;
		input >> z;
		Point P4(x, y, z); //Creates fourth point
		cout << "Point: " << P4 << "initialized!" << endl;
        out << "Point: " << P4 << "initialized!" << endl;

		Tetrahedron Te(P1, P2, P3, P4); //Creates the tetrahedron
		cout << "Tetrahedron Te initialized with: " << P1 << "" << P2 << "" << P3 << "" << P4 << endl << endl;

		int choice;
		input >> choice;

		switch (choice) { //Reads which from the file which operation should be done and does it if it is valid
		case 1: //Checks if the tetrahedron is regular
			{
				if (Te.check_if_Regular())
                {
                    cout << "Tetrahedron is regular" << endl; //Prints the result
   					out << "Tetrahedron is regular" << endl; //Prints the result
                }
				else
                {
 					cout << "Tetrahedron is irregular" << endl; //Prints the result
					out << "Tetrahedron is irregular" << endl; //Prints the result
                }
				break;
			}
			case 2: //Checks if the tetrahedron is orthogonal
			{
				if (Te.check_if_Ortogonal())
                {
 					cout << "Tetrahedron is orthogonal" << endl; //Prints the result
					out << "Tetrahedron is orthogonal" << endl; //Prints the result
                }
				else
                {
 					cout << "Tetrahedron is not orthogonal" << endl; //Prints the result
					out << "Tetrahedron is not orthogonal" << endl; //Prints the result
                }
				break;
			}
			case 3: //Finds the surface of the tetrahedron
			{
				cout << "Surface of the Tetrahedron is " << Te.get_Surface() << endl; //Prints the result
				out << "Surface of the Tetrahedron is " << Te.get_Surface() << endl; //Prints the result
				break;
			}
			case 4: //Finds the volume of the tetrahedron
			{
				cout << "Volume of the Tetrahedron is " << Te.get_Volume() << endl; //Prints the result
				out << "Volume of the Tetrahedron is " << Te.get_Volume() << endl; //Prints the result
				break;
			}
			case 5: //Checks if a point lays on the surface of the tetrahedron and is possessed inside of it
			{
				cout << "\tEnter values for x, y, z to define the Point P" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (Te.operator<(P))
                {
                    cout << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed INSIDE it" << endl;  //Prints the result
                    out << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed INSIDE it" << endl;  //Prints the result
                }
				else{
					cout << "Point " << P << " DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					cout << "  is not possessed INSIDE it or both" << endl;  //Prints the result
					out << "Point " << P << " DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					out << "  is not possessed INSIDE it or both" << endl;
				}
				break;
			}
			case 6: //Checks if a point lays on the surface of the tetrahedron and is possessed outside of it
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (Te.operator>(P))
                {
                    cout << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed OUTSIDE of it" << endl;  //Prints the result
					out << "Point " << P << " LAYS on the surface of the Tetrahedron and is possessed OUTSIDE of it" << endl;  //Prints the result
                }
				else{
					cout << "Point " << P << "DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					cout << " is not possessed OUTSIDE of it or both" << endl;  //Prints the result
					out << "Point " << P << "DOES NOT LAY on the surface of the Tetrahedron or" << endl;
					out << " is not possessed OUTSIDE of it or both" << endl;
				}
				break;
			}
			case 7: //Checks if a point lays on the segment of the tetrahedron
			{
				cout << "\tEnter values for x, y, z to define the starting Point P2" << endl;
				cout << "\tx = ";
				cin >> x;
				cout << "\ty = ";
				cin >> y;
				cout << "\tz = ";
				cin >> z;
				Point P(x, y, z); //Creates the point
				cout << "Point: " << P << "initialized!" << endl;
				out << "Point: " << P << "initialized!" << endl;

				if (Te.operator==(P))
                {
                    cout << "Point " << P << " lays on a Segment of the Tetrahedron" << endl; //Prints the result
                    out << "Point " << P << " lays on a Segment of the Tetrahedron" << endl; //Prints the result
                }
				else
                {
 					cout << "Point " << P << " does not lay on a Segment of the Triangle" << endl; //Prints the result
                    out << "Point " << P << " does not lay on a Segment of the Triangle" << endl; //Prints the result
                }
				break;
			}
			default: //If input is invalid
            {
                cout << "Invalid choice!" << endl;
                out << "Invalid choice!" << endl;
                break;
            }
        }
			cout << "If you want to make another operation enter 'y' otherwise enter 'n' " << endl;
			char ch;
			cin >> ch;
			if (ch == 'y')
				other_operation = true; //Goes back to tetrahedron menu
			else{
				other_operation = false; //Breaks the while loop
				cout << "If you want to switch the object with other press 'Y' otherwise enter 'N' " << endl;
				char c;
				cin >> c;
				if (c == 'Y')
					other_object = true; //Goes back to the main menu
				else
					other_object = false; //Exit program
			}
            out.close();
            if(other_object == false)
            {
                cout << "If you want to save operations on file, press 1, otherwise press 0" << endl;
                char write;
                cin >> write;
                if(write != '1')
                {
                    remove("output.txt");
                }
        }
    }
}
