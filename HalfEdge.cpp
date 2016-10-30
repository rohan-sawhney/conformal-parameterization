#include "HalfEdge.h"
#include "Vertex.h"

double HalfEdge::cotan() const
{
    if (onBoundary) return 0.0;
    
    Eigen::Vector3d p0 = vertex->position;
    Eigen::Vector3d p1 = next->vertex->position;
    Eigen::Vector3d p2 = next->next->vertex->position;
    
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p2 - p0;
    
    return v1.dot(v2) / v1.cross(v2).norm();
}

double HalfEdge::angle() const
{
    if (onBoundary) return 0.0;
    
    Eigen::Vector3d p0 = vertex->position;
    Eigen::Vector3d p1 = next->vertex->position;
    Eigen::Vector3d p2 = next->next->vertex->position;
    
    Eigen::Vector3d v1 = p0 - p2; v1.normalize();
    Eigen::Vector3d v2 = p1 - p2; v2.normalize();
    
    return acos(fmax(-1.0, fmin(1.0, v1.dot(v2))));
}
