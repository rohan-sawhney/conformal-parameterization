#include "Mesh.h"
#include "MeshIO.h"
#include <Eigen/SparseCholesky>
#define EPSILON 1e-6

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

void Mesh::buildConformalEnergy(Eigen::SparseMatrix<std::complex<double>>& E) const
{
    std::vector<Eigen::Triplet<std::complex<double>>> ETriplets;
    
    // build dirichlet energy
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double sumCoefficients = 0.0;
        do {
            // (cotA + cotB) / 2
            double coefficient = 0.5 * (he->cotan() + he->flip->cotan());
            sumCoefficients += coefficient;
            
            ETriplets.push_back(Eigen::Triplet<std::complex<double>>(v->index,
                                                                     he->flip->vertex->index,
                                                                     coefficient));
            
            he = he->flip->next;
        } while (he != v->he);
        
        ETriplets.push_back(Eigen::Triplet<std::complex<double>>(v->index, v->index,
                                                                 -sumCoefficients + EPSILON));
    }

    // subtract area term from dirichlet energy
    std::complex<double> i(0, 0.5);
    for (std::vector<HalfEdgeIter>::const_iterator it = boundaries.begin(); it != boundaries.end(); it++) {
        HalfEdgeCIter he = *it;
        
        do {
            int id1 = he->vertex->index;
            int id2 = he->flip->vertex->index;
            
            ETriplets.push_back(Eigen::Triplet<std::complex<double>>(id1, id2, i));
            ETriplets.push_back(Eigen::Triplet<std::complex<double>>(id2, id1, -i));
            
            he = he->next;
        } while (he != *it);
    }

    E.setFromTriplets(ETriplets.begin(), ETriplets.end());
}

void Mesh::buildMassMatrix(Eigen::SparseMatrix<std::complex<double>>& M) const
{
    std::vector<Eigen::Triplet<std::complex<double>>> MTriplets;
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        MTriplets.push_back(Eigen::Triplet<std::complex<double>>(v->index, v->index, v->dualArea()));
    }
    
    M.setFromTriplets(MTriplets.begin(), MTriplets.end());
}

double residual(const Eigen::SparseMatrix<std::complex<double>>& A,
                const Eigen::SparseMatrix<std::complex<double>>& M,
                const Eigen::VectorXcd& x)
{
    Eigen::VectorXcd b = A*x;
    std::complex<double> lambda = x.adjoint().dot(b) / x.adjoint().dot(M*x);
    return (b - lambda*M*x).norm() / x.norm();
}

void solveInversePowerMethod(const Eigen::SparseMatrix<std::complex<double>>& A,
                             const Eigen::SparseMatrix<std::complex<double>>& M,
                             const Eigen::VectorXcd& id, Eigen::VectorXcd& x)
{
    // uses inverse power method to compute smallest eigenvalue
    // subject to constraints <x, 1> = 0 and |x| = 1
    
    // prefactor
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<std::complex<double>>> solver(A);
    const Eigen::VectorXcd Mid = M * id;
    
    while (residual(A, M, x) > EPSILON) {
        // backsolve
        x = solver.solve(M*x);
        
        // center
        x -= std::conj(x.adjoint().dot(Mid))*id;
        
        // normalize
        x /= sqrt(std::norm(x.adjoint().dot(M*x)));
    }
}

void Mesh::parameterize()
{
    int v = (int)vertices.size();
    
    // build conformal energy
    Eigen::SparseMatrix<std::complex<double>> E(v, v);
    buildConformalEnergy(E);
    
    // build mass matrix
    Eigen::SparseMatrix<std::complex<double>> M(v, v);
    buildMassMatrix(M);

    // find eigenvector corresponding to smallest eigenvalue
    Eigen::VectorXcd id = Eigen::VectorXcd::Ones(v); id /= sqrt(std::norm(id.adjoint().dot(M*id)));
    Eigen::VectorXcd z = Eigen::VectorXcd::Random(v);
    solveInversePowerMethod(E, M, id, z);

    // set uv coords
    Eigen::Vector2d avg = Eigen::Vector2d::Zero();
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->uv(0) = z(v->index).real();
        v->uv(1) = z(v->index).imag();
        avg += v->uv;
    }
    avg /= v;
    
    // compute areas
    double area = 0.0;
    double uvArea = 0.0;
    for (FaceCIter f = faces.begin(); f != faces.end(); f++) {
        area += f->area();
        uvArea += f->uvArea();
    }
    
    // normalize uvs
    double scale = sqrt(area / uvArea);
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->uv = scale*(v->uv - avg);
    }
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
    }
    
    // determine radius
    double rMax = 0;
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
