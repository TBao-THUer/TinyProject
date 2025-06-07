#include    <bits/stdc++.h>
#include    "Vector.h"
#include    "Matrix.h"
#include    "LinearSystem.h"
#include    "GeneralLinearSystem.h"
using namespace std;

// Function: Load dataset from filename
void LoadDataset(const string& filename, Matrix& A, Vector& b, int& numSamples) {
    ifstream file(filename);
    string line;
    vector<vector<double>> features;
    vector<double> targets;

    // Read data
    while (getline(file, line)) {
        stringstream ss(line);
        string token;
        vector<double> row;
        int col = 0;

        while (getline(ss, token, ',')) {
            if (col >= 2 && col <= 7) { // MYCT to CHMAX (features)
                row.push_back(stod(token));
            } else if (col == 8) { // PRP (target)
                targets.push_back(stod(token));
            }
            col++;
        }

        if (row.size() == 6) {
            features.push_back(row);
        }
    }

    // Declare a matrix contains features and a vector contains targets
    numSamples = features.size();
    A = Matrix(numSamples, 6);
    b = Vector(numSamples);

    // Assign features and targets
    for(int i=1; i<=numSamples; i++) {
        for(int j=1; j<=6; j++) {
            A(i, j) = features[i - 1][j - 1];
        }
        b[i] = targets[i - 1];
    }

}

// Compute the root mean square error between predictions and targets
double ComputeRMSE(const Vector& preds, const Vector& targets) {
    int size = preds.size();
    double error = 0.0;
    for(int i=1; i<=size; i++) {
        double diff = preds[i] - targets[i];
        error += diff*diff;
    }
    return sqrt(error/size);
}

int main() {
    Matrix A_full;
    Vector b_full;
    int totalSamples;

    // Step 1: Load dataset
    LoadDataset("machine.data", A_full, b_full, totalSamples);

    // Step 2: Randomly split into 80% training and 20% validation
    vector<int> indices(totalSamples);
    for(int i=0; i<totalSamples; i++) {
        indices[i] = i+1;
    }

    mt19937 g(time(0)); // Random seed each time
    shuffle(indices.begin(), indices.end(), g);

    int train_size = static_cast<int>(0.8 * totalSamples);
    int eval_size = totalSamples - train_size;

    Matrix A_train(train_size, 6);
    Vector b_train(train_size);
    Matrix A_eval(eval_size, 6);
    Vector b_eval(eval_size);

    for(int i=1; i<=train_size; i++) {
        for(int j=1; j<=6; j++) {
            A_train(i, j) = A_full(indices[i-1], j);
        }
        b_train[i] = b_full[indices[i-1]];
    }

    for(int i=1; i<=eval_size; i++) {
        for(int j=1; j<=6; j++) {
            A_eval(i, j) = A_full(indices[i+train_size-1], j);
        }
        b_eval[i] = b_full[indices[i+train_size-1]];
    }

    // Step 3: Solve for weights using least squares
    GeneralLinearSystem system(&A_train, &b_train);
    Vector x = system.SolveMoorePenrose();

    // Step 4: Predict on test set
    Vector y_pred(eval_size);
    for(int i=1; i<=eval_size; i++) {
        double pred = 0.0;
        for(int j=1; j<=6; j++) {
            pred += A_eval(i, j) * x[j];
        }
        y_pred[i] = pred;
    }

    // Step 5: Evaluate
    double rmse = ComputeRMSE(y_pred, b_eval);
    cout << "Model coefficients (x1 to x6):\n";
    for(int i=1; i<=6; i++) {
        cout << x[i] << (i < 6 ? " " : "\n"); // -0.00904481 0.0188597 0.00230201 0.790052 2.60229 -0.200239 
    }

    cout << "Root mean square error: " << rmse << "\n";  // 151.432
    return 0;
}