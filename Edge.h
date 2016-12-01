#ifndef EDGE_H
#define EDGE_H

#include "Types.h"

class Edge {
public:
    // one of the two half edges associated with this edge
    HalfEdgeIter he;
    
    // id between 0 and |E|-1
    int index;
    
    // checks if this edge lies on boundary
    bool isBoundary() const;
    
    // returns length
    double length() const;
  
    // returns cotan weight
    double cotanWeigth() const;
    
    // flip edge
    void flip();
};

#endif
