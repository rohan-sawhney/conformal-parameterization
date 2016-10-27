#ifndef SCP_H
#define SCP_H

#include "Parameterization.h"

class Scp: public Parameterization {
public:
    // constructor
    Scp(Mesh& mesh0);
    
    // parameterize
    void parameterize() override;
    
private:
    // EC = ED - A
    void buildConformalEnergy(Eigen::SparseMatrix<std::complex<double>>& E) const;
    
    // set uvs
    void setUvs(const Eigen::VectorXcd& z);
};

#endif 
