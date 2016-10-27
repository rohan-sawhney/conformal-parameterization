#ifndef QC_ERROR_H
#define QC_ERROR_H

#include "Types.h"

class QuasiConformalError {
public:
    // computes the quasi-conformal error in a triangle with
    // vertices (p1, p2, p3) and texture coordinates (q1, q2, q3)
    static double compute(std::vector<Eigen::Vector3d> p, std::vector<Eigen::Vector3d> q);
    
    // standard color map for quasi-conformal error
    static Eigen::Vector3d color(double qc);
    
private:
    // returns hsv values
    static Eigen::Vector3d hsv(double h, double s, double v);
};

#endif
