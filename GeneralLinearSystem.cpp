#include "GeneralLinearSystem.h"
#include <bits/stdc++.h>
using namespace std;

// Parametered Constructor
GeneralLinearSystem::GeneralLinearSystem(Matrix* pA, Vector* pb) : mpA(pA), mpb(pb) {
    if (!pA || !pb) throw invalid_argument("Null pointer in GeneralLinearSystem Constructor.");
    mRows = pA->getNumRows();
    mCols = pA->getNumCols();

    if ((mRows != pb->size()) && (mCols != pb->size())) {
        throw invalid_argument("Incompatible matrix/vector dimensions.");
    }
}

// Destructor
GeneralLinearSystem::~GeneralLinearSystem() {}

// Solve using Moore-Penrose pseudoinverse
Vector GeneralLinearSystem::SolveMoorePenrose() {
    Matrix A_T = mpA->transpose();
    Matrix pseudoInv;

    if (mRows >= mCols) {
        // Overdetermined: A+ = (A^T A)^-1 A^T
        Matrix ATA = A_T * (*mpA);      // matrix * matrix
        Matrix ATA_inv = ATA.inverse();
        pseudoInv = ATA_inv * A_T;
    } else {
        // Underdetermined: A+ = A^T (A A^T)^-1
        Matrix AAT = (*mpA) * A_T;
        Matrix AAT_inv = AAT.inverse();
        pseudoInv = A_T * AAT_inv;
    }
    return pseudoInv * (*mpb);
}

// Solve using Tikhonov regularization
Vector GeneralLinearSystem::SolveTikhonov(double lambda) {
    Matrix A_T = mpA->transpose();
    Matrix regularized, invTerm;

    if (mRows >= mCols) {
        // Overdetermined: (A^T A + λI)^-1 A^T b
        Matrix ATA = A_T * (*mpA);
        Matrix I = Matrix::Identity(mCols);
        regularized = ATA + (I * lambda);
        invTerm = regularized.inverse();
        return invTerm * A_T * (*mpb);
    } else {
        // Underdetermined: A^T (A A^T + λI)^-1 b
        Matrix AAT = (*mpA) * A_T;
        Matrix I = Matrix::Identity(mRows);
        regularized = AAT + (I * lambda);
        invTerm = regularized.inverse();
        return A_T * (invTerm * (*mpb));
    }
}
