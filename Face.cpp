#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"

bool Face::isBoundary() const
{
    return he->onBoundary;
}

Eigen::Vector3d Face::normal() const
{
    Eigen::Vector3d a = he->vertex->position;
    Eigen::Vector3d b = he->next->vertex->position;
    Eigen::Vector3d c = he->next->next->vertex->position;
    
    Eigen::Vector3d v1 = b - a;
    Eigen::Vector3d v2 = c - a;
    
    return v1.cross(v2);
}

double Face::area() const
{
    if (isBoundary()) {
        return 0;
    }
    
    return 0.5 * normal().norm();
}