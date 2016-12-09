#ifndef SOLVER_H
#define SOLVER_H

#include "Types.h"

struct MeshHandle {
    // typedefs
    typedef std::function<void(double&, const Eigen::VectorXd&)> ComputeObjective;
    typedef std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> ComputeGradient;
    typedef std::function<void(Eigen::SparseMatrix<double>&, const Eigen::VectorXd&)> ComputeHessian;
    
    // constructor
    MeshHandle() {}
    
    // member variables
    ComputeObjective computeObjective;
    ComputeGradient computeGradient;
    ComputeHessian computeHessian;
};

class Solver {
public:
    // constructor
    Solver(int n0);

    // gradient descent 
    void gradientDescent();

    // coordinate descent
    void coordinateDescent();

    // newton
    void newton();

    // lbfgs
    void lbfgs();

    // member variables
    MeshHandle *handle;
    Eigen::VectorXd x;

private:
    // member variable
    int n;
};

#endif
