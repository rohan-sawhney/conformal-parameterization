#include "Lscm.h"
#include <eigen/SparseCholesky>

Lscm::Lscm(Mesh& mesh0):
Parameterization(mesh0)
{
    
}

void Lscm::pinVertices()
{
    pinnedVertices.resize(2);
    pinnedPositions = {Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 0)};
    
     // pin diameter vertices on longest boundary loop
    double max = 0.0;
    for (std::vector<HalfEdgeIter>::const_iterator it1 = mesh.boundaries.begin();
                                                   it1 != mesh.boundaries.end();
                                                   it1++) {
        HalfEdgeCIter he1 = *it1;
        do {
            const int& vIdx1(he1->vertex->index);
            const Eigen::Vector3d& p1(he1->vertex->position);
            
            for (std::vector<HalfEdgeIter>::const_iterator it2 = mesh.boundaries.begin();
                                                           it2 != mesh.boundaries.end();
                                                           it2++) {
                HalfEdgeCIter he2 = *it2;
                do {
                    const int& vIdx2(he2->vertex->index);
                    const Eigen::Vector3d& p2(he2->vertex->position);
                    double l = (p2-p1).squaredNorm();
                    
                    if (l > max) {
                        max = l;
                        pinnedVertices[0] = 2*vIdx1;
                        pinnedVertices[1] = 2*vIdx2;
                    }
                    
                    he2 = he2->next;
                } while (he2 != *it2);
            }
            
            he1 = he1->next;
        } while (he1 != *it1);
    }
}

void computeLocalBasisCoordinates(std::vector<Eigen::Vector2d>& coords,
                                  const std::vector<Eigen::Vector3d>& positions)
{
    Eigen::Vector3d x = positions[1] - positions[0];
    Eigen::Vector3d y = positions[2] - positions[0];
    
    // compute orthonormal basis
    Eigen::Vector3d xhat = x; xhat.normalize();
    Eigen::Vector3d zhat = xhat.cross(y); zhat.normalize();
    Eigen::Vector3d yhat = zhat.cross(xhat); yhat.normalize();
    
    // compute coordinates in local basis
    coords.push_back(Eigen::Vector2d(0, 0));
    coords.push_back(Eigen::Vector2d(x.norm(), 0));
    coords.push_back(Eigen::Vector2d(y.dot(xhat), y.dot(yhat)));
}

bool Lscm::isPinnedVertex(const int& vIndex, int& shift, Eigen::Vector2d& pinnedPosition) const
{
    shift = 0;
    for (size_t i = 0; i < pinnedVertices.size(); i++) {
        if (pinnedVertices[i] == vIndex) {
            pinnedPosition = pinnedPositions[i];
            return true;
        }
        
        if (pinnedVertices[i] < vIndex) shift += 2;
    }
    
    return false;
}

void Lscm::buildMassMatrix(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& b) const
{
    std::vector<Eigen::Triplet<double>> MTriplets;
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fIndex = 2*f->index;
            
            // get indices and positions of adjacent vertices
            std::vector<int> indices;
            std::vector<Eigen::Vector3d> positions;
            HalfEdgeCIter h = f->he;
            do {
                indices.push_back(2*h->vertex->index);
                positions.push_back(h->vertex->position);
                h = h->next;
                
            } while (h != f->he);
            
            // compute weights
            std::vector<Eigen::Vector2d> coords;
            computeLocalBasisCoordinates(coords, positions);
            std::vector<Eigen::Vector2d> ws = {coords[2] - coords[1],
                                               coords[0] - coords[2],
                                               coords[1] - coords[0]};
            
            // set M entries
            for (int i = 0; i < 3; i++) {
                int shift;
                Eigen::Vector2d pinnedPosition;
                const Eigen::Vector2d& w(ws[i]);
                
                if (isPinnedVertex(indices[i], shift, pinnedPosition)) {
                    b(fIndex) -= (w.x()*pinnedPosition.x() - w.y()*pinnedPosition.y());
                    b(fIndex+1) -= (w.y()*pinnedPosition.x() + w.x()*pinnedPosition.y());
                    
                } else {
                    int vIndex = indices[i] - shift;
                    
                    // set real components
                    MTriplets.push_back(Eigen::Triplet<double>(fIndex, vIndex, w.x()));
                    MTriplets.push_back(Eigen::Triplet<double>(fIndex+1, vIndex+1, w.x()));
                    
                    if (i != 2) {
                        // set imaginary components
                        MTriplets.push_back(Eigen::Triplet<double>(fIndex, vIndex+1, -w.y()));
                        MTriplets.push_back(Eigen::Triplet<double>(fIndex+1, vIndex, w.y()));
                    }
                }
            }
        }
    }
    
    M.setFromTriplets(MTriplets.begin(), MTriplets.end());
}

void solveLeastSquares(const Eigen::SparseMatrix<double>& A,
                       const Eigen::VectorXd& b,
                       Eigen::VectorXd& x)
{
    // initialize
    Eigen::SparseMatrix<double> At = A.transpose();
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(At * A);
    Eigen::VectorXd Atb = At * b;
    
    // backsolve
    x = solver.solve(Atb);
}

void Lscm::setUvs(const Eigen::VectorXd& x)
{
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int vIndex = 2*v->index;
        
        int shift;
        Eigen::Vector2d pinnedPosition;
        if (isPinnedVertex(vIndex, shift, pinnedPosition)) {
            v->uv(0) = pinnedPosition.x();
            v->uv(1) = pinnedPosition.y();
            
        } else {
            vIndex -= shift;
            v->uv(0) = x(vIndex);
            v->uv(1) = x(vIndex+1);
        }
    }
    
    normalize();
}

void Lscm::parameterize()
{
    int f = 2*(int)(mesh.faces.size() - mesh.boundaries.size());
    int v = 2*((int)mesh.vertices.size() - 2);
    
    // pin vertices
    pinVertices();
    
    // build mass matrix and handle boundary conditions
    Eigen::SparseMatrix<double> M(f, v);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(f);
    buildMassMatrix(M, b);
    
    // solve least squares
    Eigen::VectorXd x(v);
    solveLeastSquares(M, b, x);
    
    // set uv coords
    setUvs(x);
}
