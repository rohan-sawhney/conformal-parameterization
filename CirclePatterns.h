#ifndef CIRCLE_PATTERNS_H
#define CIRCLE_PATTERNS_H

#include "Parameterization.h"
#include "MosekSolver.h"
#include <stack>

class CirclePatterns: public Parameterization {
public:
    // constructor
    CirclePatterns(Mesh& mesh0);
    
    // parameterize
    void parameterize() override;
    
protected:
    // sets angle based constraints, bounds and objective function
    void setupAngleOptProblem();
        
    // sets thetas
    void setThetas();
    
    // compute angles
    bool computeAngles();
    
    // sets radii based constraints and bounds
    void setupRadiiOptProblem();
    
    // computes radii energy, its gradient and hessian, and their sparsity
    void computeEnergy(double& energy, const double *rho);
    void computeGradient(double *gradient, const double *rho);
    void computeHessian(double *hessian, const double *rho);
    void buildGradientSparsity(int *idx);
    void buildHessianSparsity(int *idxi, int *idxj);

    // sets radii
    void setRadii();
    
    // compute radii
    bool computeRadii();
    
    // computes angles and edge lengths
    void computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths);
    
    // determines position of unfixed face vertex
    void performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir, Eigen::VectorXd& lengths,
                           std::unordered_map<int, bool>& visited, std::stack<EdgeCIter>& stack);
    
    // sets uvs
    void setUVs();
    
    // member variables
    Eigen::VectorXd angles;
    Eigen::VectorXd thetas;
    Eigen::VectorXd radii;
    Eigen::VectorXi eIntIndices;
    int imaginaryHe;
    MosekSolver solver;
};

#endif 
