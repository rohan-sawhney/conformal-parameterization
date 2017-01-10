#include "Cetm.h"

Cetm::Cetm(Mesh& mesh0, int optScheme0):
Parameterization(mesh0),
thetas(mesh.vertices.size()),
lengths(mesh.edges.size()),
angles(mesh.halfEdges.size()),
OptScheme(optScheme0)
{
    
}

void Cetm::setTargetThetas()
{
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (!v->isBoundary()) thetas[index[v->index]] = 2*M_PI;
    }
}

double angle(double lij, double ljk, double lki)
{
    return acos(fmax(-1.0, fmin(1.0, (lij*lij + lki*lki - ljk*ljk) / (2*lij*lki))));
}

double lambda(double l)
{
    return 2*log(l);
}

double lobachevsky(double a)
{
    return Cl2(2*a)/2;
}

double cot(double a)
{
    // handle degenerate case by clamping (see Section 3.1 of CETM paper)
    if (a == 0.0 || a == M_PI) return 0.0;
    return tan(M_PI_2 - a);
}

void Cetm::computeEnergy(double& energy, const Eigen::VectorXd& u)
{
    energy = 0.0;
    
    // sum over faces
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            
            // copy the three u values for triangle f into au
            int i = 0;
            std::vector<double> au(3);
            HalfEdgeCIter he = f->he;
            do {
                if (he->vertex->isBoundary()) au[i] = 0.0;
                else au[i] = u[index[he->vertex->index]];
                i++;
                
                he = he->next;
            } while (he != f->he);
            
            // copy the three edge lengths values for triangle f into al
            i = 0;
            std::vector<double> al(3);
            do {
                al[i] = he->edge->length()*exp(0.5*(au[i] + au[(i+1)%3]));
                i++;
                
                he = he->next;
            } while (he != f->he);
            
            // add f(~λij, ~λjk, ~λk) - π(ui + uj + uk)/2 term to energy
            i = 0;
            do {
                angles[he->index] = angle(al[(i+2)%3], al[i], al[(i+1)%3]);
                energy += 0.5*angles[he->index]*lambda(al[i]) + lobachevsky(angles[he->index]) - M_PI_2*au[i];
                i++;
                
                he = he->next;
            } while (he != f->he);
        }
    }
    
    // sum over vertices
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->isBoundary()) continue;
        int vIdx = index[v->index];
        energy += 0.5*thetas[vIdx]*u[vIdx];
    }
}

void Cetm::computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& u)
{
    // loop over vertices
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->isBoundary()) continue;
        
        // compute angle sum
        double angleSum = 0.0;
        HalfEdgeCIter he = v->he;
        do {
            if (!he->onBoundary) angleSum += angles[he->next->index];
            
            he = he->flip->next;
        } while (he != v->he);
        
        int vIdx = index[v->index];
        gradient[vIdx] = 0.5*(thetas[vIdx] - angleSum);
    }
}

void Cetm::computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u)
{
    std::vector<Eigen::Triplet<double>> HTriplets;
    
    // build dirichlet energy
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->isBoundary()) continue;
        
        int vIdx = index[v->index];
        HalfEdgeCIter he = v->he;
        double sumW = 0.0;
        do {
            // (cotA + cotB) / 4
            double cotAlpha = !he->onBoundary ? cot(angles[he->index]) : 0.0;
            double cotBeta  = !he->flip->onBoundary ? cot(angles[he->flip->index]) : 0.0;
            double w = (cotAlpha + cotBeta)/4.0;
            sumW += w;
            
            if (!he->flip->vertex->isBoundary()) {
                HTriplets.push_back(Eigen::Triplet<double>(vIdx, index[he->flip->vertex->index], -w));
            }
            
            he = he->flip->next;
        } while (he != v->he);
        
        HTriplets.push_back(Eigen::Triplet<double>(vIdx, vIdx, sumW + 1e-8));
    }
    
    hessian.setFromTriplets(HTriplets.begin(), HTriplets.end());
}

void Cetm::setEdgeLengthsAndAngles()
{
    // set edge lengths
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        double u1 = he->vertex->isBoundary() ? 0.0 : solver.x[index[he->vertex->index]];
        double u2 = he->flip->vertex->isBoundary() ? 0.0 : solver.x[index[he->flip->vertex->index]];
        
        lengths[e->index] = e->length()*exp(0.5*(u1 + u2));
    }
    
    // set angles
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            
            // copy the three edge lengths values for triangle f into al
            int i = 0;
            std::vector<double> al(3);
            HalfEdgeCIter he = f->he;
            do {
                al[i] = lengths[he->edge->index];
                i++;
                
                he = he->next;
            } while (he != f->he);
            
            // compute angles
            i = 0;
            do {
                angles[he->index] = angle(al[(i+2)%3], al[i], al[(i+1)%3]);
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
    if (OptScheme == GRAD_DESCENT) solver.gradientDescent();
    else if (OptScheme == NEWTON) solver.newton();
    else if (OptScheme == TRUST_REGION) solver.trustRegion();
    else solver.lbfgs();
    
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
    // Set vertex indices
    int vIdx = 0;
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (!v->isBoundary()) index[v->index] = vIdx++;
    }
    solver.n = vIdx;
    
    // set target thetas
    setTargetThetas();
    
    // compute edge lengths
    if (!computeScaleFactors()) {
        std::cout << "Unable to compute edge lengths" << std::endl;
        return;
    }

    // set uvs
    setUVs();
}
