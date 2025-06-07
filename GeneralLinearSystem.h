#ifndef GENERAL_LINEARSYSTEM_H
#define GENERAL_LINEARSYSTEM_H

#include    "Matrix.h"
#include    "Vector.h"

class GeneralLinearSystem {
private:
    Matrix* mpA;
    Vector* mpb;
    int mRows, mCols;

public:
    // Constructor
    GeneralLinearSystem(Matrix* pA, Vector* pb);
    // Destructor
    ~GeneralLinearSystem();

    Vector SolveMoorePenrose();                 // Least-squares/min-norm solution
    Vector SolveTikhonov(double lambda=1e-6);   // Regularized solution
};

#endif // GENERAL_LINEARSYSTEM_H