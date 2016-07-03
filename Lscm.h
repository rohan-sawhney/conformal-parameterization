#ifndef LSCM_H
#define LSCM_H

#include "Parameterizer.h"

class Lscm : public Parameterizer {
public:
    // constructor
    Lscm(Mesh& mesh0);
    
    // parameterize
    void parameterize();
    
private:
    // finds diameter vertices
    void findDiameterVertices();
    
    // checks if vertex is pinned
    bool isPinnedVertex(const int& vIndex, int& shift, Eigen::Vector2d& pinnedPosition) const;
    
    // builds mass matrix
    void buildMassMatrix(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& b) const;
    
    // set uvs
    void setUvs(const Eigen::VectorXd& x);
    
    // member variables
    std::vector<Eigen::Vector2d> pinnedPositions;
    std::vector<int> pinnedVertices;
};

#endif
