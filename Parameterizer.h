#ifndef PARAMETERIZER_H
#define PARAMETERIZER_H

#include "Mesh.h"
#include <Eigen/SparseCore>

class Parameterizer {
public:
    // constructor
    Parameterizer(Mesh& mesh0);
    
    // destructor
    virtual ~Parameterizer() {}
    
    // flatten
    virtual void parameterize() = 0;
    
protected:
    // normalize
    void normalize();
    
    // member variable
    Mesh& mesh;
};

#endif 
