#include "Parameterizer.h"

Parameterizer::Parameterizer(Mesh& mesh0):
mesh(mesh0)
{
    
}

void Parameterizer::normalize()
{
    // compute areas
    double area = 0.0;
    double uvArea = 0.0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        area += f->area();
        uvArea += f->uvArea();
    }
    
    // compute center
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        center += v->uv;
    }
    center /= (double)mesh.vertices.size();
    
    // shfit and scale
    double scale = sqrt(area / uvArea);
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv = scale*(v->uv - center);
    }
}