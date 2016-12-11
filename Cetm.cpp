#include "Cetm.h"
#define cot(x) (tan(M_PI_2 - x))

Cetm::Cetm(Mesh& mesh0):
Parameterization(mesh0),
thetas(mesh.vertices.size()),
lengths(mesh.edges.size()),
angles(mesh.halfEdges.size()),
solver((int)mesh.vertices.size())
{
    
}

void Cetm::setDefaultThetas()
{
    // set target thetas to a disk
    int imaginaryHe = 0;
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (e->isBoundary()) imaginaryHe++;
    }
    
    double angleSum = M_PI * (imaginaryHe - 2) / imaginaryHe;;
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->isBoundary()) thetas[v->index] = angleSum;
        else thetas[v->index] = 2*M_PI;
    }
}

double angle(double lij, double ljk, double lki)
{
    return acos((lij*lij + lki*lki - ljk*ljk) / (2*lij*lki));
}

bool validAngle(double a)
{
    return a > 0 && a < M_PI;
}

double lambda(double l)
{
    return 2*log(l);
}

double lobachevsky(double a)
{
    return Cl2(2*a)/2;
}

void Cetm::computeEnergy(double& energy, const Eigen::VectorXd& u)
{
    energy = 0.0;
    
    // sum over faces
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            
            // add adjacent u values
            int i = 0;
            std::vector<double> au;
            HalfEdgeCIter he = f->he;
            do {
                au.push_back(u[he->vertex->index]);
                energy -= M_PI_2*au[i];
                i++;
                
                he = he->next;
            } while (he != f->he);
            
            // add adjacent edge lengths
            i = 0;
            std::vector<double> al;
            do {
                al.push_back(he->edge->length() * exp(0.5*(au[i] + au[(i+1)%3])));
                i++;
                
                he = he->next;
            } while (he != f->he);
            
            // add angles
            i = 0; int l = -1;
            do {
                // check if triangle inequality is not satisfied
                int j = (i+1)%3, k = (i+2)%3;
                if (al[i] > al[k] + al[j]) {
                    l = i;
                    break;
                }
                
                angles[he->index] = angle(al[k], al[i], al[j]);
                i++;
                
                he = he->next;
            } while (he != f->he);
            
            // add f(λij', λjk', λk') term to energy
            i = 0;
            he = f->he;
            do {
                if (l != -1) angles[he->index] = l == i ? M_PI - EPSILON : EPSILON/2.0;
                energy += 0.5*angles[he->index]*lambda(al[i]) + lobachevsky(angles[he->index]);
                i++;
                
                he = he->next;
            } while (he != f->he);
        }
    }
    
    // sum over vertices
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        energy += 0.5*thetas[v->index]*u[v->index];
    }
}

void Cetm::computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& u)
{
    // loop over vertices
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        // compute angle sum
        double angleSum = 0.0;
        HalfEdgeCIter he = v->he;
        do {
            if (!he->onBoundary) angleSum += angles[he->next->index];
            
            he = he->flip->next;
        } while (he != v->he);
        
        gradient[v->index] = 0.5*(thetas[v->index] - angleSum);
    }
}

void Cetm::computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u)
{
    std::vector<Eigen::Triplet<double>> HTriplets;
    
    // build dirichlet energy
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double sumW = 0.0;
        do {
            // (cotA + cotB) / 4
            double w = (!he->onBoundary && validAngle(angles[he->index])) ?
                        cot(angles[he->index]) : 0.0;
            w += (!he->flip->onBoundary && validAngle(angles[he->flip->index])) ?
                  cot(angles[he->flip->index]) : 0.0;
            w /= 4.0;
            sumW += w;
            
            HTriplets.push_back(Eigen::Triplet<double>(v->index, he->flip->vertex->index, -w));
            
            he = he->flip->next;
        } while (he != v->he);
        
        HTriplets.push_back(Eigen::Triplet<double>(v->index, v->index, sumW + 1e-8));
    }
    
    hessian.setFromTriplets(HTriplets.begin(), HTriplets.end());
}

void Cetm::setEdgeLengthsAndAngles()
{
    // set edge lengths
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        lengths[e->index] = e->length() * exp(0.5*(solver.x[he->vertex->index] +
                                                   solver.x[he->flip->vertex->index]));
    }
    
    // set angles
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            
            // add adjacent edge lengths
            std::vector<double> al;
            HalfEdgeCIter he = f->he;
            do {
                al.push_back(lengths[he->edge->index]);
                
                he = he->next;
            } while (he != f->he);
            
            int i = 0; int l = -1;
            do {
                // check if triangle inequality is not satisfied
                int j = (i+1)%3, k = (i+2)%3;
                if (al[i] > al[k] + al[j]) {
                    l = i;
                    break;
                }
                
                angles[he->index] = angle(al[k], al[i], al[j]);
                i++;
                
                he = he->next;
            } while (he != f->he);
            
            i = 0;
            he = f->he;
            do {
                if (l != -1) angles[he->index] = l == i ? M_PI - EPSILON : EPSILON/2.0;
                i++;
                
                he = he->next;
            } while (he != f->he);
        }
    }
}

bool Cetm::computeScaleFactors()
{
    MeshHandle handle;
    handle.computeEnergy = std::bind(&Cetm::computeEnergy, this, _1, _2);
    handle.computeGradient = std::bind(&Cetm::computeGradient, this, _1, _2);
    handle.computeHessian = std::bind(&Cetm::computeHessian, this, _1, _2);
    
    solver.handle = &handle;
    solver.newton();
    
    // set edge lengths and angles
    setEdgeLengthsAndAngles();
    
    return true;
}

void Cetm::performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir,
                             std::unordered_map<int, bool>& visited, std::stack<EdgeCIter>& stack)
{
    if (!he->onBoundary) {
        int fIdx = he->face->index;
        if (visited.find(fIdx) == visited.end()) {
            HalfEdgeCIter next = he->next;
            HalfEdgeCIter prev = he->next->next;
            
            // compute new uv position
            double angle = angles[next->index];
            Eigen::Vector2d newDir = {cos(angle)*dir[0] - sin(angle)*dir[1],
                                      sin(angle)*dir[0] + cos(angle)*dir[1]};
            
            prev->vertex->uv = he->vertex->uv + newDir*lengths[prev->edge->index];
            
            // mark face as visited
            visited[fIdx] = true;
            
            // push edges onto stack
            stack.push(next->edge);
            stack.push(prev->edge);
        }
    }
}

void Cetm::setUVs()
{
    // push any edge
    std::stack<EdgeCIter> stack;
    EdgeCIter e = mesh.edges.begin();
    stack.push(e);
    e->he->vertex->uv = Eigen::Vector2d::Zero();
    e->he->next->vertex->uv = Eigen::Vector2d(lengths[e->index], 0);
    
    // perform layout
    std::unordered_map<int, bool> visited;
    while (!stack.empty()) {
        EdgeCIter e = stack.top();
        stack.pop();
        
        HalfEdgeCIter h1 = e->he;
        HalfEdgeCIter h2 = h1->flip;
        
        // compute edge vector
        Eigen::Vector2d dir = h2->vertex->uv - h1->vertex->uv;
        dir.normalize();
        
        performFaceLayout(h1, dir, visited, stack);
        performFaceLayout(h2, -dir, visited, stack);
    }
    
    normalize();
}

void Cetm::parameterize()
{
    // set default thetas
    setDefaultThetas();
    
    // compute edge lengths
    if (!computeScaleFactors()) {
        std::cout << "Unable to compute edge lengths" << std::endl;
        return;
    }

    // set uvs
    setUVs();
    
    // TODOs
    // 1) implement natural boundary conditions
}
