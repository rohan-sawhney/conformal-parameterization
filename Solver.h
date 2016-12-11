#ifndef SOLVER_H
#define SOLVER_H

#include "Types.h"
using namespace std::placeholders;

struct MeshHandle {
    // typedefs
    typedef std::function<void(double&, const Eigen::VectorXd&)> ComputeEnergy;
    typedef std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> ComputeGradient;
    typedef std::function<void(Eigen::SparseMatrix<double>&, const Eigen::VectorXd&)> ComputeHessian;
    
    // constructor
    MeshHandle() {}
    
    // member variables
    ComputeEnergy computeEnergy;
    ComputeGradient computeGradient;
    ComputeHessian computeHessian;
};

class Solver {
public:
    // constructor
    Solver(int n0);
    
    // gradient descent
    void gradientDescent();
    
    // newton
    void newton();
    
    // lbfgs
    void lbfgs(int m = 10);
    
    // member variables
    MeshHandle *handle;
    Eigen::VectorXd x;
    
private:
    // member variable
    int n;
};

#endif
