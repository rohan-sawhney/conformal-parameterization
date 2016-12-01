#include "Subdivision.h"

Eigen::Vector3d Subdivision::updateVertexPosition(const VertexCIter& v) const
{
    Eigen::Vector3d coordinate;
    if (v->isBoundary()) {
        coordinate = 3.0*v->position/4.0;
        HalfEdgeCIter h = v->he;
        do {
            if (h->onBoundary) {
                HalfEdgeCIter prev = h;
                do {
                    prev = prev->next;
                } while (prev->next != h);
                
                coordinate += (h->flip->vertex->position + prev->vertex->position)/8.0;
                break;
            }
            
            h = h->flip->next;
        } while (h != v->he);
        
    } else {
        int n = 0;
        Eigen::Vector3d sum = Eigen::Vector3d::Zero();
        HalfEdgeCIter h = v->he;
        do {
            sum += h->flip->vertex->position;
            n++;
            
            h = h->flip->next;
        } while (h != v->he);
        
        double beta = n == 3 ? 3.0/16.0 : 3.0/(8.0*n);
        coordinate = (1 - n*beta)*v->position + beta*sum;
    }
    
    return coordinate;
}

Eigen::Vector3d Subdivision::createNewCoordinate(const HalfEdgeCIter& h) const
{
    Eigen::Vector3d coordinate;
    VertexCIter v1 = h->vertex;
    VertexCIter v2 = h->flip->vertex;
    
    if (h->flip->onBoundary) {
        coordinate = (v1->position + v2->position)/2.0;
        
    } else {
        VertexCIter v3 = h->next->next->vertex;
        VertexCIter v4 = h->flip->next->next->vertex;
        
        coordinate = (v1->position*3.0 + v2->position*3.0 +
                      v3->position + v4->position)/8.0;
    }
    
    return coordinate;
}

void Subdivision::createNewFaces(MeshData& data, const std::vector<int>& indices) const
{
    std::vector<Index> faceIndices;
    faceIndices.push_back(Index(indices[0],-1,-1));
    faceIndices.push_back(Index(indices[1],-1,-1));
    faceIndices.push_back(Index(indices[5],-1,-1));
    data.indices.push_back(faceIndices);
    
    faceIndices.clear();
    faceIndices.push_back(Index(indices[5],-1,-1));
    faceIndices.push_back(Index(indices[1],-1,-1));
    faceIndices.push_back(Index(indices[3],-1,-1));
    data.indices.push_back(faceIndices);
    
    faceIndices.clear();
    faceIndices.push_back(Index(indices[3],-1,-1));
    faceIndices.push_back(Index(indices[1],-1,-1));
    faceIndices.push_back(Index(indices[2],-1,-1));
    data.indices.push_back(faceIndices);
    
    faceIndices.clear();
    faceIndices.push_back(Index(indices[5],-1,-1));
    faceIndices.push_back(Index(indices[3],-1,-1));
    faceIndices.push_back(Index(indices[4],-1,-1));
    data.indices.push_back(faceIndices);
}

void Subdivision::subdivide(Mesh& mesh) const
{
    MeshData data;
    
    // insert updated position of existing vertices
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        data.positions.push_back(updateVertexPosition(v));
    }
    
    // create new coordinates and faces
    std::unordered_map<int, bool> visited;
    std::unordered_map<int, int> edgeVertexMap;
    data.indices.reserve(4*mesh.faces.size());
    for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            std::vector<int> indices;
            HalfEdgeCIter h = f->he;
            do {
                indices.push_back(h->vertex->index);
                
                int eIdx = h->edge->index;
                if (!visited[eIdx]) {
                    edgeVertexMap[eIdx] = (int)data.positions.size();
                    data.positions.push_back(createNewCoordinate(h));
                    
                    visited[eIdx] = true;
                }
                indices.push_back(edgeVertexMap[eIdx]);
                
                h = h->next;
            } while (h != f->he);
            
            createNewFaces(data, indices);
        }
    }
    
    // build mesh
    MeshIO::buildMesh(data, mesh);
}

