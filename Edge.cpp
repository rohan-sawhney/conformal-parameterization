#include "Edge.h"
#include "HalfEdge.h"
#include "Vertex.h"
#include "Face.h"

double Edge::length() const
{
    Eigen::Vector3d a = he->vertex->position;
    Eigen::Vector3d b = he->flip->vertex->position;
    
    return (b-a).norm();
}

bool Edge::isBoundary() const
{
    return he->onBoundary || he->flip->onBoundary;
}

double Edge::cotanWeigth() const
{
    return 0.5 * (he->cotan() + he->flip->cotan());
}

void Edge::flip()
{
    HalfEdgeIter heNext = he->next;
    HalfEdgeIter heNextNext = heNext->next;
    
    HalfEdgeIter flip = he->flip;
    HalfEdgeIter flipNext = flip->next;
    HalfEdgeIter flipNextNext = flipNext->next;
    
    VertexIter v1 = he->vertex;
    VertexIter v2 = heNext->vertex;
    VertexIter v3 = heNextNext->vertex;
    VertexIter v4 = flipNextNext->vertex;
    
    FaceIter f = he->face;
    FaceIter fFlip = flip->face;
    
    // set halfEdge vertex
    he->vertex = v4;
    flip->vertex = v3;
    
    // set next halfEdge
    he->next = heNextNext;
    heNextNext->next = flipNext;
    flipNext->next = he;
    
    flip->next = flipNextNext;
    flipNextNext->next = heNext;
    heNext->next = flip;
    
    // set halfEdge face
    heNext->face = fFlip;
    flipNext->face = f;
    
    // set face halEdge
    f->he = he;
    fFlip->he = flip;
    
    // set vertex halfEdge
    if (v1->he == he) v1->he = flipNext;
    if (v2->he == flip) v2->he = heNext;
}
