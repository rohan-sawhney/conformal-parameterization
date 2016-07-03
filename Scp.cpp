#include "Scp.h"
#include <Eigen/SparseCholesky>
#define MAX_ITER 32

Scp::Scp(Mesh& mesh0):
Parameterizer(mesh0)
{
    
}

void Scp::buildConformalEnergy(Eigen::SparseMatrix<std::complex<double>>& E) const
{
    std::vector<Eigen::Triplet<std::complex<double>>> ETriplets;
    
    // build dirichlet energy
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double sumCoefficients = 0.0;
        do {
            // (cotA + cotB) / 4
            double coefficient = 0.25 * (he->cotan() + he->flip->cotan());
            sumCoefficients += coefficient;
            
            ETriplets.push_back(Eigen::Triplet<std::complex<double>>(v->index,
                                                                     he->flip->vertex->index,
                                                                     -coefficient));
            
            he = he->flip->next;
        } while (he != v->he);
        
        ETriplets.push_back(Eigen::Triplet<std::complex<double>>(v->index, v->index,
                                                                 sumCoefficients));
    }
    
    // subtract area term from dirichlet energy
    std::complex<double> i(0, 0.25);
    for (std::vector<HalfEdgeIter>::const_iterator it = mesh.boundaries.begin();
                                                   it != mesh.boundaries.end();
                                                   it++) {
        HalfEdgeCIter he = *it;
        do {
            int id1 = he->vertex->index;
            int id2 = he->flip->vertex->index;
            
            ETriplets.push_back(Eigen::Triplet<std::complex<double>>(id1, id2, -i));
            ETriplets.push_back(Eigen::Triplet<std::complex<double>>(id2, id1, i));
            
            he = he->next;
        } while (he != *it);
    }
    
    E.setFromTriplets(ETriplets.begin(), ETriplets.end());
}

void Scp::buildMassMatrix(Eigen::SparseMatrix<std::complex<double>>& M) const
{
    std::vector<Eigen::Triplet<std::complex<double>>> MTriplets;
    
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        MTriplets.push_back(Eigen::Triplet<std::complex<double>>(v->index, v->index, v->dualArea()));
    }
    
    M.setFromTriplets(MTriplets.begin(), MTriplets.end());
}

double residual(const Eigen::SparseMatrix<std::complex<double>>& A,
                const Eigen::VectorXcd& x)
{
    Eigen::VectorXcd b = A*x;
    std::complex<double> lambda = x.adjoint().dot(b) / x.adjoint().dot(x);
    return (b - lambda*x).norm() / x.norm();
}

void solveInversePowerMethod(const Eigen::SparseMatrix<std::complex<double>>& A,
                             Eigen::VectorXcd& x)
{
    // uses inverse power method to compute smallest eigenvalue
    // subject to constraints <x, 1> = 0 and |x| = 1
    
    // prefactor
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<std::complex<double>>> solver(A);
    
    for (int i = 0; i < MAX_ITER; i++) {
        // backsolve
        x = solver.solve(x);
        
        // center
        std::complex<double> mean = x.sum() / (double)x.size();
        for (int i = 0; i < x.size(); i++) x(i) -= mean;
        
        // normalize
        x.normalize();
    }
    
    std::cout << "residual: " << residual(A, x) << std::endl;
}

void Scp::setUvs(const Eigen::VectorXcd& z)
{
    // set uv coords
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        v->uv(0) = z(v->index).real();
        v->uv(1) = z(v->index).imag();
    }
    
    normalize();
}

void Scp::parameterize()
{
    int v = (int)mesh.vertices.size();
    
    // build conformal energy
    Eigen::SparseMatrix<std::complex<double>> E(v, v);
    buildConformalEnergy(E);
    
    // build mass matrix
    Eigen::SparseMatrix<std::complex<double>> M(v, v);
    buildMassMatrix(M);
    
    E += std::complex<double>(EPSILON)*M;
    
    // find eigenvector corresponding to smallest eigenvalue
    Eigen::VectorXcd z = Eigen::VectorXcd::Random(v);
    solveInversePowerMethod(E, z);
    
    // set uv coords
    setUvs(z);
}
