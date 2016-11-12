#ifndef PARAMETERIZER_H
#define PARAMETERIZER_H

#include "Mesh.h"

class Parameterization {
public:
    // constructor
    Parameterization(Mesh& mesh0);
    
    // destructor
    virtual ~Parameterization() {}
    
    // flatten
    virtual void parameterize() = 0;
    
    // computes quasi conformal error
    double computeQcError();
    
protected:
    // normalize
    void normalize();
    
    // member variable
    Mesh& mesh;
};

#endif 
