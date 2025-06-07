#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include    "Vector.h"
#include    "Matrix.h"

class LinearSystem {
protected:
    int mSize;
    Matrix *mpA;
    Vector *mpb;

private:
    // Prevent use of automatically generated copy constructor
    LinearSystem(const LinearSystem& other) = delete;
    // Prevent use of automatically generated assignment operator
    LinearSystem& operator=(const LinearSystem& other) = delete;

public:
    // Constructor
    LinearSystem(Matrix* pA, Vector* pb);

    // Destructor
    virtual ~LinearSystem();

    // Virtual solve method to be overriden by derived classes
    virtual Vector Solve();

    // Function GaussianElimination
    void GaussianElimination(Vector& solution);
};

class PosSymLinSystem : public LinearSystem {
public:
    // Constructor
    PosSymLinSystem(Matrix* pA, Vector *pb);

    // Overriden Solve method to use conjugate gradient
    Vector Solve() override;

private:
    // Check if matrix is symmetric
    bool IsSymmetric() const;
};

#endif  // LINEARSYSTEM_H