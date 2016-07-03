#ifndef SCP_H
#define SCP_H

#include "Parameterizer.h"

class Scp : public Parameterizer {
public:
    // constructor
    Scp(Mesh& mesh0);
    
    // parameterize
    void parameterize();
    
private:
    // EC = ED - A
    void buildConformalEnergy(Eigen::SparseMatrix<std::complex<double>>& E) const;
    
    // builds mass matrix
    void buildMassMatrix(Eigen::SparseMatrix<std::complex<double>>& M) const;
    
    // set uvs
    void setUvs(const Eigen::VectorXcd& z);
};

#endif 
