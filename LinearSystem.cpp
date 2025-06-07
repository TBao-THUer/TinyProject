#include    "LinearSystem.h"
#include    <bits/stdc++.h>
using namespace std;

// LinearSystem constructor
LinearSystem::LinearSystem(Matrix* pA, Vector* pb) : mpA(pA), mpb(pb) {
    if (pA == nullptr || pb == nullptr) {
        throw invalid_argument("Null pointer passed to LinearSystem constructor.");
    }

    if (pA->getNumCols() != pA->getNumRows()) {
        throw invalid_argument("Matrix must be square.");
    }

    if (pA->getNumRows() != pb->size()) {
        throw invalid_argument("Matrix and vector sizes are incompatible.");
    }
    mSize = pA->getNumRows();
}

// Destructor
LinearSystem::~LinearSystem() {}

// Solving using Gaussian elimination with pivoting
Vector LinearSystem::Solve() {
    Vector solution(mSize);
    GaussianElimination(solution);
    return solution;
}

// Gaussian elimination with partial pivoting
void LinearSystem::GaussianElimination(Vector& solution) {
    Matrix A(*mpA);     // Copy constructor 
    Vector b(*mpb);     // Copy constructor

    for(int k=1; k<=mSize; k++) {
        // Partial pivoting
        int maxRow = k;
        double maxVal = abs(A(k, k));
        for(int i=k+1; i<=mSize; i++) {
            if (abs(A(i, k)) > maxVal) { 
                maxVal = abs(A(i, k));
                maxRow = i;
            }
        }

        // Swap rows if necessary
        if (maxRow != k) {
            for(int j=k; j<=mSize; j++) {
                swap(A(k, j), A(maxRow, j));
            }
            swap(b[k], b[maxRow]);
        }

        // Elimination
        for(int i=k+1; i<=mSize; i++) {
            double factor = A(i, k) / A(k, k);
            for(int j=k; j<=mSize; j++) {
                A(i, j) -= factor * A(k, j);
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    for(int i=mSize; i>0; i--) {
        solution[i] = b[i];
        for(int j=i+1; j<=mSize; j++) {
            solution[i] -= A(i, j) * solution[j];
        }
        solution[i] /= A(i, i);
    }
}

// PosSymLinSystem constructor
PosSymLinSystem::PosSymLinSystem(Matrix* pA, Vector* pb) : LinearSystem(pA, pb) {
    if (!IsSymmetric()) {
        throw invalid_argument("Matrix is not symmetric.");
    }
}

// Solve using conjugate gradient method
Vector PosSymLinSystem::Solve() {
    const double tolerance = 1e-10;
    const int maxIterations = 1000;

    Vector x(mSize); // initial guess: zero
    Vector r = (*mpb) - (*mpA) * x;
    Vector p = r;

    double rsold = r * r;

    for (int i = 0; i < maxIterations; ++i) {
        Vector Ap = (*mpA) * p;
        double alpha = rsold / (p * Ap);

        x = x + alpha * p;
        r = r - alpha * Ap;

        double rsnew = r * r;
        if (sqrt(rsnew) < tolerance) {
            break;
        }

        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    }

    return x;
}

// Check for symmetry
bool PosSymLinSystem::IsSymmetric() const {
    const double symmetryTolerance = 1e-10;

    for (int i = 1; i <= mSize; ++i) {
        for (int j = 1; j < i; ++j) {
            if (std::abs((*mpA)(i, j) - (*mpA)(j, i)) > symmetryTolerance) {
                return false;
            }
        }
    }

    return true;
}