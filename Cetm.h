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
    
    // computes energy, gradient and hessian
    void computeEnergy(double& energy, const Eigen::VectorXd& u);
    void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& u);
    void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u);
    
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
    Solver solver;
};

#endif
