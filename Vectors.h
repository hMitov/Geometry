#ifndef VECTORS_H_INCLUDED
#define VECTORS_H_INCLUDED

#include <exception>
#include <cstring>
/*
class Element
{
public:
    virtual void fun() = 0;
};
*/

class Point//: public Element
{
protected:
    float x;
    float y;
    float z;
public:
    Point(); // default constructor
    Point(float x, float y, float z); // constructor
    ~Point(); // destructor
    Point(const Point& rhs); // copy constructor
    Point& operator=(const Point& rhs); // operator=

    bool operator==(Point& rhs); // checks if two points coincide
    friend std::ostream& operator<<(std::ostream& out, const Point& rhs);

    friend class Vector;
    friend class Line;
    friend class Segment;
    friend class Triangle;
    friend class Tetrahedron;
};

class VectorLengthException: public std::exception
{
public:
    const char* what()const throw()
    {
        return "Error! An attempt to find the NULL vector length!";
    }
}obj;


class EqualPointException: public std::exception
{
public:
    const char* what()const throw()
    {
        return "Error! An attempt to initialize Triangle with three equal points!" ;
    }
    const char* whatt()const throw()
    {
        return "Error! An attempt to initialize Tetrahedron with two or more equal points!";
    }

}obj1;

class Vector: public Point
{
protected:
    float x, y, z;
public:
    friend class Line;
    Vector(); // default constructor
    Vector(float x, float y, float z); // constructor
    Vector(float x, float y, float z, Point* p1, Point* p2); // constructor
    ~Vector(); // destructor
    Vector(const Vector& rhs); // copy constructor
    Vector& operator=(const Vector& rhs); // operator=


    float vector_length(); // finds length of vector
    Vector vector_direction(); // finds vector direction
    bool zero_vector_check(); // checks if zero vector
    bool parallel_vectors(Vector& v); // checks if vectors are parallel
    bool perpendicular_vectors(Vector& v); // checks if vectors are perpendicular

    Vector operator+(const Vector& v); // finds sum of two vectors
    Vector operator-(const Vector& v); // finds subtraction of two vectors
    float operator*(const Vector& v); // finds vector-multiply of two vectors
    Vector operator*(int num); // finds multiply of vector and number
    Vector operator^(const Vector& v); // finds scalar multiply of two vectors
    float operator()( Vector& v,  Vector& w); // finds mixed multiply
/*
    friend class Segment;
    friend class Triangle;*/

    friend std::ostream& operator<<(std::ostream& out, const Vector& rhs);
};

class Line: public Vector
{
protected:
    Point starting_point;
    Vector line_direction;
public:
    Line(); // default constructor
    Line(float x, float y, float z, Point& p1, Vector& v); // constructor
    Line(float x, float y, float z, Point& p1, Point& p2); // constructor
    ~Line(); // destructor
    Line(const Line& rhs); // copy constructor
    Line& operator=(const Line& rhs); // operator=

    Vector find_direction(); // finds direction of line
    Vector normal_vector(); // finds the normal vector
    float Angel_between_two_lines(Line& rhs); // finds the angle between two lines
    bool operator+(Point& p); // checks if Point& p lays on the line
    bool operator||(Line& rhs); // checks if two lines are parallel
    bool operator==(Line& rhs); // checks if Line& rhs coincides with implicit line
    bool operator&&(Line& l); //checks if Line& l comes across with implicit line
    bool operator!=(Line& l); //checks if Line& l is crossed with implicit line
    bool operator|(Line& l); //checks if two lines are perpendicular
};

class Segment: public Line
{
protected:
    Point start_point;
    Point end_point;
public:
    Segment(); // default constructor
    Segment(Point& p1, Vector& v, Point& p2); // constructor
    ~Segment(); // destructor
    Segment(const Segment& rhs); // copy constructor
    Segment& operator=(const Segment& rhs); // operator=

    float find_length_of_Segment(); // finds the length of the segment
    Point find_Mid_point_of_Segment(); // finds the middle point of the segment
    bool operator==(Point& p); // checks if a point lays on the segment
};


class Triangle: public Point
{
protected:
    Point A;
    Point B;
    Point C;
public:
    Triangle(); // default constructor
    Triangle(float x, float y, float z, Point& p1, Point& p2, Point& p3); // constructor
    ~Triangle(); // destructor
    Triangle(const Triangle& rhs); // copy constructor
    Triangle& operator=(const Triangle& rhs); // operator=

    void define_Type(); // finds the type of the triangle(ravnostranen...ostrougulen)
    double find_Area();
    double find_Perimeter();
    Point find_Medicenter();
    float does_point_lays_on_(Point& p); // checks if a point lays on the triangle

    bool operator<(Point& p); // checks if Point& p lays on the surface of the triangle and is inner
    bool operator>(Point& p); // checks if Point& p lays on the surface of the triangle and is outer
    bool operator==(Point& p); // checks if Point& p lays at least on one line of the triangle
};



class Tetrahedron: public Point
{
protected:
    Point p1;
    Point p2;
    Point p3;
    Point p4;
public:
    Tetrahedron(); // default constructor
    Tetrahedron(float x, float y, float z, Point& p1, Point& p2, Point& p3, Point& p4); // constructor
    ~Tetrahedron(); // destructor
    Tetrahedron(const Tetrahedron& rhs); // copy constructor
    Tetrahedron& operator=(const Tetrahedron& rhs); // operator=

    bool check_if_Regular(); // checks if all lines are equal
    bool check_if_Ortogonal();
    float get_Surface(); // surface area
    float get_Volume(); // obem

    bool operator<(Point& p); // checks if Point& p lays on the surface of the tetrahedron and is inner
    bool operator>(Point& p); // checks if Point& p lays on the surface of the tetrahedron and is outer
    bool operator==(Point& p); // checks if Point& p lays at least on one line of the tetrahedron
};



class Menu {
    char choose;
public:
    void mainMenu(); // this function calls all below object-Menu functions

    static bool readData(char choose); // in this function we enter either we want to read from console or from file.
                                       // Our preference is stored in "choose" variable

    void pointMenu(); // initializes object of class Point and contains all operations available for that class

    void vectorMenu(); // initializes object of class Vector and contains all operations available for that class

    void lineMenu(); // initializes object of class Line and contains all operations available for that class

    void segmentMenu(); // initializes object of class Segment and contains all operations available for that class

    void triangleMenu(); // initializes object of class Triangle and contains all operations available for that class

    void tetrahedronMenu(); // initializes object of class Tetrahedron and contains all operations available for that class

    static bool other_object; // static variable that stores either we want program to ask us for another object
                              // after ending the current operation or not

    static bool other_operation;  // static variable that stores either we want program to ask us for another operation
                                  // after ending the current operation or not. If not asks us for another object.
};

bool Menu:: other_object = true;
bool Menu:: other_operation = true;


#endif // VECTORS_H_INCLUDED
