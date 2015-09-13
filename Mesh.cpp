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
            // (cotA + cotB) / 4
            double coefficient = 0.25 * (he->cotan() + he->flip->cotan());
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
    std::complex<double> i(0, 0.25);
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

double residual(const Eigen::SparseMatrix<std::complex<double>>& A, const Eigen::VectorXcd& x)
{
    Eigen::VectorXcd b = A*x;
    std::complex<double> lambda = x.dot(b);
    return (b - lambda*x).norm();
}

void solveInversePowerMethod(const Eigen::SparseMatrix<std::complex<double>>& A, Eigen::VectorXcd& x)
{
    // uses inverse power method to compute smallest eigenvalue
    // subject to constraints <x, 1> = 0 and |x| = 1
    
    // prefactor
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<std::complex<double>>> solver(A);
    int s = (int)x.size();
    
    while (residual(A, x) > EPSILON) {
        // backsolve
        x = solver.solve(x);
        
        // subtract mean
        std::complex<double> mean = x.sum() / (double)s;
        for (int i = 0; i < s; i++) x(i) -= mean;
        
        // normalize
        x.normalize();
    }
}

void Mesh::parameterize()
{
    int v = (int)vertices.size();
    
    // build conformal energy
    Eigen::SparseMatrix<std::complex<double>> E(v, v);
    buildConformalEnergy(E);

    // find eigenvector corresponding to smallest eigenvalue
    Eigen::VectorXcd z = Eigen::VectorXcd::Random(v);
    solveInversePowerMethod(E, z);

    // set uv coords
    double max = 0.0;
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->uv(0) = z(v->index).real();
        v->uv(1) = z(v->index).imag();
        max = std::max(max, v->uv.norm());
    }
    
    // normalize uvs
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->uv /= max;
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
