#include "Parameterization.h"
#include "QcError.h"

Parameterization::Parameterization(Mesh& mesh0):
mesh(mesh0)
{
    
}

double uvArea(FaceCIter f)
{
    if (f->isBoundary()) {
        return 0;
    }
    
    const Eigen::Vector2d& a(f->he->vertex->uv);
    const Eigen::Vector2d& b(f->he->next->vertex->uv);
    const Eigen::Vector2d& c(f->he->next->next->vertex->uv);
    
    const Eigen::Vector2d u = b - a;
    const Eigen::Vector2d v = c - a;
    
    return 0.5 * (u.x()*v.y() - v.x()*u.y());
}

Eigen::Vector2d uvBarycenter(FaceCIter f)
{
    if (f->isBoundary()) {
        return Eigen::Vector2d::Zero();
    }
    
    const Eigen::Vector2d& a(f->he->vertex->uv);
    const Eigen::Vector2d& b(f->he->next->vertex->uv);
    const Eigen::Vector2d& c(f->he->next->next->vertex->uv);
    
    return (a + b + c) / 3.0;
}

void Parameterization::normalize()
{
    // compute center
    double totalArea = 0;
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        double area = uvArea(f);
        center += area * uvBarycenter(f);
        totalArea += area;
    }
    center /= totalArea;
    
    // shift
    double r = 0.0;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv -= center;
        r = std::max(r, v->uv.squaredNorm());
    }
    
    // scale
    r = sqrt(r);
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv /= r;
    }
}

void Parameterization::computeQcError()
{
    for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        std::vector<Eigen::Vector3d> p, q;
        HalfEdgeCIter he = f->he;
        do {
            p.push_back(he->vertex->position);
            q.push_back(Eigen::Vector3d(he->vertex->uv.x(), he->vertex->uv.y(), 0));
            
            he = he->next;
        } while (he != f->he);
        
        f->qcError = QuasiConformalError::color(QuasiConformalError::compute(p, q));
    }
}
