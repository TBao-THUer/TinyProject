#include    "Vector.h"
#include    <bits/stdc++.h>
using namespace std;

// Default constructor
Vector::Vector() : mSize(0), mData(nullptr) {}

// Parametered constructor
Vector::Vector(int size): mSize(size)
{
    if (size <= 0) throw invalid_argument("Size must be positive.");
    mData = new double[size];
    for(int i=0; i<size; i++)
    {
        mData[i] = 0.0;
    }
}

// Copy constructor
Vector::Vector(const Vector& other) : mSize(other.mSize) 
{
    mData = new double[mSize];
    for(int i=0; i<mSize; i++)
    {
        mData[i] = other.mData[i];
    }
}

// Destructor
Vector::~Vector()
{
    delete[] mData;
}

// Assignment operator
Vector& Vector::operator=(const Vector& other)
{
    if (this == &other) return *this; // self-assignment check

    delete[] mData;
    mSize = other.mSize;
    mData = new double[mSize];
    for(int i=0; i<mSize; i++)
    {
        mData[i] = other.mData[i];
    }
    return *this;
}

// Unary Minus
Vector Vector::operator-() const
{
    Vector ans(mSize);
    for(int i=0; i<mSize; i++)
    {
        ans.mData[i] = -mData[i];
    }
    return ans;
}

// Addition operator
Vector Vector::operator+(const Vector& other) const
{
    if (mSize != other.mSize) throw length_error("Size mismatch.");
    Vector ans(mSize);
    for(int i=0; i<mSize; i++)
    {
        ans.mData[i] = mData[i] + other.mData[i];
    }
    return ans;
}

// Substraction
Vector Vector::operator-(const Vector& other) const
{
    if (mSize != other.mSize) throw length_error("Size mismatch.");
    Vector ans(mSize);
    for(int i=0; i<mSize; i++)
    {
        ans.mData[i] = mData[i] - other.mData[i];
    }
    return ans;
}

// Dot product (vector * vector -> scalar)
double Vector::operator*(const Vector& other) const {
    if (mSize != other.mSize) {
        throw length_error("Dot product error: vectors must be the same size.");
    }

    double result = 0.0;
    for(int i=0; i<mSize; i++) {
        result += mData[i] * other.mData[i];
    }
    return result;
}

// Scalar Multiplication (vector * scalar)
Vector Vector::operator*(double scalar) const
{
    Vector ans(mSize);
    for(int i=0; i<mSize; i++)
    {
        ans.mData[i] = mData[i] * scalar;
    }
    return ans;
}

// Scalar Multiplication (scalar * vector)
Vector operator*(double scalar, const Vector& v)
{
    return v * scalar;
}

// Indexing Operators (1-based indexing)
double& Vector::operator[](int index)
{
    if (index < 1 || index > mSize) throw out_of_range("Index out of range (must be 1-based).");
    return mData[index-1];
}

const double& Vector::operator[](int index) const
{
    if (index < 1 || index > mSize) throw out_of_range("Index out of range (must be 1-based).");
    return mData[index - 1];
}

double& Vector::operator()(int index) 
{
    return (*this)[index];
}

const double& Vector::operator()(int index) const
{
    return (*this)[index];
}

// Get size
int Vector::size() const
{
    return mSize;
}

// Print all elements of vector
void Vector::print() const
{
    for(int i=0; i<mSize; i++)
    {
        cout << fixed << setprecision(2) << mData[i] << (i<mSize-1 ? " " : "\n");
    }
}