#include "Mesh.h"
#include "MeshIO.h"
#include "Scp.h"
#include "Lscm.h"
#include "QcError.h"

Mesh::Mesh()
{
    
}

Mesh::Mesh(const Mesh& mesh)
{
    *this = mesh;
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        normalize();
    }
    
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    return false;
}

void Mesh::setQcError()
{
    for (FaceIter f = faces.begin(); f != faces.end(); f++) {
        std::vector<Eigen::Vector3d> p, q;
        HalfEdgeCIter he = f->he;
        do {
            p.push_back(he->vertex->position);
            q.push_back(Eigen::Vector3d(he->vertex->uv.x(), he->vertex->uv.y(), 0));
            
            he = he->next;
        } while (he != f->he);
        
        f->error = QuasiConformalError::color(QuasiConformalError::compute(p, q));
    }
}

void Mesh::parameterize(const int& technique)
{
    if (technique == SCP) {
        Scp scp(*this);
        scp.parameterize();
    
    } else {
        Lscm lscm(*this);
        lscm.parameterize();
    }
    
    setQcError();
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin and determine radius
    double rMax = 0;
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
