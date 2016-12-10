#ifndef CIRCLE_PATTERNS_H
#define CIRCLE_PATTERNS_H

#include "Parameterization.h"
#include "Utils.h"
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
    
    // computes radii energy, its gradient and hessian for our solver
    void computeEnergy(double& energy, const Eigen::VectorXd& rho);
    void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho);
    void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& rho);

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
    MosekSolver::Solver mosekSolver;
    Solver solver;
};

#endif 
