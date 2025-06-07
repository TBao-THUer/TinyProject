#ifndef MATRIX_H
#define MATRIX_H

#include    <iostream>
#include    "Vector.h"
#include    <stdexcept>

class Matrix{
private:
    int mNumRows;
    int mNumCols;
    double **mData;

public:
    // Constructors and Destructor
    Matrix();                           // Default constructor
    Matrix(int numRows, int numCols);   // Parameterized constructor
    Matrix(const Matrix& other);        // Copy constructor
    ~Matrix();                          // Destructor

    // Get the mNumRows and mNumCols
    int getNumRows() const;
    int getNumCols() const;

    // Overloaded operators
    double& operator()(int i, int j);       // For assignment: A(i, j) = value;
    double operator()(int i, int j) const;  // For assignment: value = A(i, j);
    Matrix& operator=(const Matrix& other); // Assignment operator

    // Unary operators
    Matrix operator-() const;

    // Binary operations
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;    // matrix * matrix
    Vector operator*(const Vector& v) const;        // matrix * vector
    Matrix operator*(double scalar) const;          // matrix * scalar

    // Determinant, inverse, pseudo-inverse
    double determinant() const;
    Matrix inverse() const;     // Find the inverse of matrix A^-1
    Matrix transpose() const;   // Transpose matrix A^T
    Matrix pseudoInverse() const;
    static Matrix Identity(int n);  // Create an identity matrix I with the size of nxn

    // Print all elements of matrix
    void print() const;
};

#endif // MATRIX_H
