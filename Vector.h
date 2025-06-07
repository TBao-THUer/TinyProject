#ifndef VECTOR_H
#define VECTOR_H    

#include    <iostream>
#include    <stdexcept>

class Vector 
{
private:
    int mSize;
    double *mData;

public:
    // Constructor and Destructor
    Vector();                           // Default constructor
    Vector(int size);                    // Parameterized constructor
    Vector(const Vector& other);        // Copy constructor
    ~Vector();                          // Destructor

    // Assignment Operator
    Vector& operator=(const Vector& other);

    // Unary Minus
    Vector operator-() const;

    // Binary Operators
    Vector operator+(const Vector& other) const; // Addition
    Vector operator-(const Vector& other) const; // Substraction

    // Dot product (vector * vector -> scalar)
    double operator*(const Vector& other) const;

    // Scalar Multiplication
    Vector operator*(double scalar) const;
    friend Vector operator*(double scalar, const Vector& v); // scalar * vector

    // Indexing Operators (1-based indexing)
    double& operator[](int index);   // square brackets
    const double& operator[](int index) const;

    double& operator()(int index);   // round brackets
    const double& operator()(int index) const;

    // Get the size of vector
    int size() const;

    // Print all elements of vector
    void print() const;
};

#endif  // VECTOR_H