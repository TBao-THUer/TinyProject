#include    "Matrix.h"
#include    <bits/stdc++.h>
using namespace std;

// Default constructor
Matrix::Matrix() : mNumRows(0), mNumCols(0), mData(nullptr) {}

// Parametered constructor
Matrix::Matrix(int numRows, int numCols) : mNumRows(numRows), mNumCols(numCols) {
    if (numRows <= 0 || numCols <= 0) throw invalid_argument("NumRows and numCols must be positive.");
    mData = new double* [numRows];
    for(int r=0; r<mNumRows; r++){
        mData[r] = new double [numCols];
        for(int c=0; c<mNumCols; c++){
            mData[r][c] = 0;
        }
    }
}

// Copy Constructor
Matrix::Matrix(const Matrix& other) : mNumRows(other.mNumRows), mNumCols(other.mNumCols){
    mData = new double* [mNumRows];
    for(int r=0; r<mNumRows; r++){
        mData[r] = new double [mNumCols];
        for(int c=0; c<mNumCols; c++){
            mData[r][c] = other.mData[r][c];
        }
    }
}

// Destructor
Matrix::~Matrix() {
    for(int i=0; i<mNumRows; i++) {
        delete[] mData[i];
    }
    delete[] mData;
}

// Get the numRows
int Matrix::getNumRows() const {
    return mNumRows;
}

// Get the numCols
int Matrix::getNumCols() const {
    return mNumCols;
}

// Round brackets (get the value of row ith and col jth 1-based)
double& Matrix::operator()(int i, int j) {
    if (i <= 0 || i > mNumRows || j <= 0 || j > mNumCols) throw length_error("The index is out of range.");
    return mData[i - 1][j - 1];
}

double Matrix::operator()(int i, int j) const {
    if (i <= 0 || i > mNumRows || j <= 0 || j > mNumCols) throw length_error("The index is out of range.");
    return mData[i - 1][j - 1];
}
Matrix& Matrix::operator=(const Matrix& other) {
    if (this == &other) return *this;

    // Free existing memory
    for(int i=0; i<mNumRows; i++) delete[] mData[i];
    delete[] mData;

    mNumRows = other.mNumRows;
    mNumCols = other.mNumCols;

    mData = new double* [mNumRows];
    for(int r=0; r<mNumRows; r++) {
        mData[r] = new double [mNumCols];
        for(int c=0; c<mNumCols; c++) {
            mData[r][c] = other.mData[r][c];
        }
    }
    return *this;
}

// Unary operator
Matrix Matrix::operator-() const {
    Matrix ans(mNumRows, mNumCols);
    for(int r=0; r<mNumRows; r++) {
        for(int c=0; c<mNumCols; c++) {
            ans.mData[r][c] = -mData[r-1][c-1];
        }
    }
    return ans;
}

// Addition
Matrix Matrix::operator+(const Matrix& other) const {
    if (mNumRows != other.mNumRows || mNumCols != other.mNumCols) throw length_error("Size mismatch.");
    Matrix ans(mNumRows, mNumCols);
    for(int r=0; r<mNumRows; r++) {
        for(int c=0; c<mNumCols; c++) {
            ans.mData[r][c] = mData[r][c] + other.mData[r][c];
        }
    }
    return ans;
}

// Substraction
Matrix Matrix::operator-(const Matrix& other) const {
    if (mNumRows != other.mNumRows || mNumCols != other.mNumCols) throw length_error("Size mismatch.");
    Matrix ans(mNumRows, mNumCols);
    for(int r=0; r<mNumRows; r++) {
        for(int c=0; c<mNumCols; c++) {
            ans.mData[r][c] = mData[r][c] - other.mData[r][c];
        }
    }
    return ans;
}

// Multiplycation
Matrix Matrix::operator*(const Matrix& other) const {
    if (mNumCols != other.mNumRows) throw length_error("Size mismatch.");
    Matrix ans(mNumRows, other.mNumCols);
    for(int i=0; i<mNumRows; i++) {
        for(int j=0; j<other.mNumCols; j++) {
            for(int k=0; k<mNumCols; k++) {
                ans.mData[i][j] += mData[i][k] * other.mData[k][j];
            }
        }
    }
    return ans;
}

// Matrix * vector
Vector Matrix::operator*(const Vector& v) const {
    if (mNumCols != v.size()) throw length_error("The size of matrix and vector mismatch.");
    Vector ans(mNumRows);
    for(int r=0; r<mNumRows; r++) {
        for(int c=0; c<mNumCols; c++) { 
            ans[r+1] += mData[r][c] * v[c+1];
        }
    }
    return ans;
}

// Matrix * scalar
Matrix Matrix::operator*(double scalar) const {
    Matrix ans(mNumRows, mNumCols);
    for(int r=0; r<mNumRows; r++) {
        for(int c=0; c<mNumCols; c++) {
            ans.mData[r][c] = mData[r][c] * scalar;
        }
    }
    return ans;
}

// Scalar Multiplication (scalar * matrix)
Matrix operator*(double scalar, const Matrix& m)
{
    return m * scalar;
}

// Get minor
Matrix getMinor(const Matrix& mat, int row, int col) {
    Matrix minor(mat.getNumRows()-1, mat.getNumCols() - 1);
    int m_row = 1;
    for (int i = 1; i <= mat.getNumRows(); ++i) {
        if (i == row) continue;
        int m_col = 1;
        for (int j = 1; j <= mat.getNumCols(); ++j) {
            if (j == col) continue;
            minor(m_row, m_col) = mat(i, j);
            ++m_col;
        }
        ++m_row;
    }
    return minor;
}

// Get the determinant of the given matrix
double Matrix::determinant() const {
    if (mNumRows != mNumCols) throw length_error("The matrix is not square.");
    if (mNumRows == 1) {
        return mData[0][0];
    }

    if (mNumRows == 2) { 
        return mData[0][0] * mData[1][1] - mData[0][1] * mData[1][0];
    }

    double det = 0.0;
    for(int c=1; c<=mNumCols; c++) {
        Matrix minor = getMinor(*this, 1, c);
        double cofactor = ((c % 2 == 1) ? 1 : -1) * (*this)(1, c);
        det += cofactor * minor.determinant();
    }
    return det;
}

// Get the inverse of matrix
Matrix Matrix::inverse() const {
    if (mNumRows != mNumCols) {
        throw std::runtime_error("Matrix must be square to compute inverse.");
    }

    int n = mNumRows;
    Matrix augmented(n, 2 * n);

    // Build augmented matrix [A | I]
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            augmented(i, j) = (*this)(i, j);
        }
        for (int j = n + 1; j <= 2 * n; ++j) {
            augmented(i, j) = (j - n == i) ? 1.0 : 0.0;
        }
    }

    // Gauss-Jordan Elimination
    for (int i = 1; i <= n; ++i) {
        double pivot = augmented(i, i);
        if (std::fabs(pivot) < 1e-9) {
            throw std::runtime_error("Matrix is singular or nearly singular (zero pivot).");
        }

        for (int j = 1; j <= 2 * n; ++j) {
            augmented(i, j) /= pivot;
        }

        for (int k = 1; k <= n; ++k) {
            if (k == i) continue;
            double factor = augmented(k, i);
            for (int j = 1; j <= 2 * n; ++j) {
                augmented(k, j) -= factor * augmented(i, j);
            }
        }
    }

    // Extract inverse matrix
    Matrix inverse(n, n);
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            inverse(i, j) = augmented(i, j + n);
        }
    }

    return inverse;
}

Matrix Matrix::transpose() const {
    Matrix ans(mNumCols, mNumRows);
    for(int r=1; r<=mNumRows; r++) {
        for(int c=1; c<=mNumCols; c++) {
            ans(c, r) = (*this)(r, c);
        }
    }
    return ans;
}

Matrix Matrix::pseudoInverse() const {
    Matrix A_T = this -> transpose();   // A^T
    Matrix ATA = A_T * (*this);         // A^T * A

    if (ATA.getNumRows() != ATA.getNumCols()) { 
        throw runtime_error("A^T * A is not square - cannot compute pseudo-inverse.");
    }

    Matrix ATA_inv = ATA.inverse();     // (A^T * A)^-1
    return ATA_inv * A_T;
}

// Create the identity matrix I with the size of nxn
Matrix Matrix::Identity(int n) {
    Matrix I(n, n);
    for(int i=1; i<=n; i++) {
        I(i, i) = 1.0;
    }
    return I;
}

// Print all elements of matrix
void Matrix::print() const {
    for(int r=1; r<=mNumRows; r++) {
        for(int c=1; c<=mNumCols; c++) {
            cout << (*this)(r, c) << (c < mNumCols ? " " : "\n");
        }
    }
}