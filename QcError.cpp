#include "QcError.h"

double QuasiConformalError::compute(std::vector<Eigen::Vector3d> p, std::vector<Eigen::Vector3d> q)
{
    // Compute edge vectors
    Eigen::Vector3d u1 = p[1] - p[0];
    Eigen::Vector3d u2 = p[2] - p[0];
    
    Eigen::Vector3d v1 = q[1] - q[0];
    Eigen::Vector3d v2 = q[2] - q[0];
    
    // Compute orthonormal bases
    Eigen::Vector3d e1 = u1; e1.normalize();
    Eigen::Vector3d e2 = (u2 - u2.dot(e1)*e1); e2.normalize();
    
    Eigen::Vector3d f1 = v1; f1.normalize();
    Eigen::Vector3d f2 = (v2 - v2.dot(f1)*f1); f2.normalize();
    
    // Project onto bases
    p[0] = Eigen::Vector3d::Zero();
    p[1] = Eigen::Vector3d(u1.dot(e1), u1.dot(e2), 0);
    p[2] = Eigen::Vector3d(u2.dot(e1), u2.dot(e2), 0);
    
    q[0] = Eigen::Vector3d::Zero();
    q[1] = Eigen::Vector3d(v1.dot(f1), v1.dot(f2), 0);
    q[2] = Eigen::Vector3d(v2.dot(f1), v2.dot(f2), 0);
    
    double A = 2.0 * u1.cross(u2).norm();
    
    Eigen::Vector3d Ss = (q[0]*(p[1].y() - p[2].y()) +
                          q[1]*(p[2].y() - p[0].y()) +
                          q[2]*(p[0].y() - p[1].y())) / A;
    Eigen::Vector3d St = (q[0]*(p[2].x() - p[1].x()) +
                          q[1]*(p[0].x() - p[2].x()) +
                          q[2]*(p[1].x() - p[0].x())) / A;
    double a = Ss.dot(Ss);
    double b = Ss.dot(St);
    double c = St.dot(St);
    double Gamma = sqrt(0.5*((a+c) + sqrt((a-c)*(a-c) + 4.0*b*b)));
    double gamma = sqrt(0.5*((a+c) - sqrt((a-c)*(a-c) + 4.0*b*b)));
    
    if (Gamma < gamma) {
        std::swap(Gamma, gamma);
    }
    
    return Gamma / gamma;
}

Eigen::Vector3d QuasiConformalError::hsv(double h, double s, double v)
{
    double r = 0, g = 0, b = 0;
    
    if (s == 0) {
        r = v;
        g = v;
        b = v;
    
    } else {
        h = (h == 1 ? 0 : h) * 6;
        
        int i = (int)floor(h);
        
        double f = h - i;
        double p = v * (1 - s);
        double q = v * (1 - (s * f));
        double t = v * (1 - s * (1 - f));
        
        switch (i) {
            case 0:
                r = v;
                g = t;
                b = p;
                break;
                
            case 1:
                r = q;
                g = v;
                b = p;
                break;
                
            case 2:
                r = p;
                g = v;
                b = t;
                break;
                
            case 3:
                r = p;
                g = q;
                b = v;
                break;
                
            case 4:
                r = t;
                g = p;
                b = v;
                break;
                
            case 5:
                r = v;
                g = p;
                b = q;
                break;
                
            default:
                break;
        }
    }
    
    return Eigen::Vector3d(r, g, b);
}

Eigen::Vector3d QuasiConformalError::color(double qc)
{
    // Clamp to range [1, 1.5]
    qc = std::max(1.0, std::min(1.5, qc));
    
    // Compute color
    return hsv((2.0 - 4.0*(qc-1.0))/3.0, 0.7, 0.65);
}
