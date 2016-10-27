#include "Vertex.h"
#include "HalfEdge.h"
#include "Face.h"

std::vector<HalfEdge> isolated;

bool Vertex::isIsolated() const
{
    return he == isolated.begin();
}

bool Vertex::isBoundary() const
{
    HalfEdgeCIter h = he;
    do {
        if (h->onBoundary) return true;
        
        h = h->flip->next;
    } while (h != he);
    
    return false;
}

int Vertex::degree() const
{
    int k = 0;
    HalfEdgeCIter h = he;
    do {
        k++;
        
        h = h->flip->next;
    } while (h != he);
    
    return k;
}

Eigen::Vector3d Vertex::normal() const
{
    Eigen::Vector3d normal = Eigen::Vector3d::Zero();
    if (isIsolated()) return normal;
    
    HalfEdgeCIter h = he;
    do {
        Eigen::Vector3d e1 = h->next->vertex->position - position;
        Eigen::Vector3d e2 = h->next->next->vertex->position - position;
        
        double d = e1.dot(e2) / sqrt(e1.squaredNorm() * e2.squaredNorm());
        if (d < -1.0) d = -1.0;
        else if (d >  1.0) d = 1.0;
        double angle = acos(d);
        
        Eigen::Vector3d n = h->face->normal();
        normal += angle * n;
        
        h = h->flip->next;
    } while (h != he);
    
    if (!normal.isZero()) normal.normalize();
    return normal;
}
