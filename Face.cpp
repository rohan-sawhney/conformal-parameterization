#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"

bool Face::isBoundary() const
{
    return he->onBoundary;
}

Eigen::Vector3d Face::normal() const
{
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    return (b-a).cross(c-a);
}

double Face::area() const
{
    if (isBoundary()) {
        return 0;
    }
    
    return 0.5 * normal().norm();
}

double Face::uvArea() const
{
    if (isBoundary()) {
        return 0;
    }
    
    const Eigen::Vector2d& a(he->vertex->uv);
    const Eigen::Vector2d& b(he->next->vertex->uv);
    const Eigen::Vector2d& c(he->next->next->vertex->uv);
    
    const Eigen::Vector2d u = b - a;
    const Eigen::Vector2d v = c - a;
    
    return 0.5 * (u.x()*v.y() - v.x()*u.y());
}