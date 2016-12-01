#ifndef CETM_H
#define CETM_H

#include "Parameterization.h"
#include "Utils.h"
#include <stack>

class Cetm: public Parameterization {
public:
    // constructor
    Cetm(Mesh& mesh0);

    // parameterize
    void parameterize() override;
    
protected:
    // sets default theta values
    void setDefaultThetas();
    
    // sets scale factor based constraints and bounds
    void setupOptProblem();
    
    // computes energy, its gradient and hessian, and their sparsity
    void computeEnergy(double& energy, const double *u);
    void computeGradient(double *gradient, const double *u);
    void computeHessian(double *hessian, const double *u);
    void buildGradientSparsity(int *idx);
    void buildHessianSparsity(int *idxi, int *idxj);
    
    // sets edge lengths
    void setEdgeLengthsAndAngles();
    
    // computes scale factors
    bool computeScaleFactors();
    
    // determines position of unfixed face vertex
    void performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir, 
                           std::unordered_map<int, bool>& visited, std::stack<EdgeCIter>& stack);
    
    // sets uvs
    void setUVs();
    
    // member variables
    Eigen::VectorXd thetas;
    Eigen::VectorXd lengths;
    Eigen::VectorXd angles;
    MosekSolver solver;
};

#endif
