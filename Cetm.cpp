#include "Cetm.h"
#define cot(x) (tan(M_PI_2 - x))

Cetm::Cetm(Mesh& mesh0):
Parameterization(mesh0),
thetas(mesh.vertices.size()),
lengths(mesh.edges.size()),
angles(mesh.halfEdges.size())
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

void Cetm::setupOptProblem()
{
    // set constraint: sum of us over all vertices equals 0
    solver.bkc[0] = MSK_BK_FX;
    solver.blc[0] = 0.0;
    solver.buc[0] = 0.0;
    
    // set linear constraint vector and bounds
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int vIdx = v->index;
        
        solver.ptrb[vIdx] = vIdx;
        solver.ptre[vIdx] = vIdx+1;
        solver.asub[vIdx] = 0;
        solver.aval[vIdx] = 1.0;
        
        solver.bkx[vIdx] = MSK_BK_FR;
        solver.blx[vIdx] = -MSK_INFINITY;
        solver.bux[vIdx] = MSK_INFINITY;
        
        solver.c[vIdx] = 1.0;
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

void Cetm::computeEnergy(double& energy, const double *u)
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
                if (l != -1) angles[he->index] = l == i ? M_PI : 0.0;
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

void Cetm::computeGradient(double *gradient, const double *u)
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

void Cetm::computeHessian(double *hessian, const double *u)
{
    // add diagonal entries
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        hessian[v->index] = 0.0;
    }
    
    // add non-diagonal entries
    int index = (int)mesh.vertices.size();
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int i = v->index;
        
        HalfEdgeCIter he = v->he;
        do {
            int j = he->flip->vertex->index;
            
            // mosek requires only the upper (or lower triangular) part of the hessian
            if (i > j) {
                // compute cotan weight
                double w = (!he->onBoundary && validAngle(angles[he->index])) ? cot(angles[he->index]) : 0.0;
                w += (!he->flip->onBoundary && validAngle(angles[he->flip->index])) ? cot(angles[he->flip->index]) : 0.0;
                w /= 4.0;
                
                hessian[i] += w;
                hessian[j] += w;
                hessian[index] = -w;
                index++;
            }
            
            he = he->flip->next;
        } while (he != v->he);
    }
}

void Cetm::buildGradientSparsity(int *idx)
{
    // set indices for nonzero elements
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int vIdx = v->index;
        idx[vIdx] = vIdx;
    }
}

void Cetm::buildHessianSparsity(int *idxi, int *idxj)
{
    // set indices for nonzero diagonal elements
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int i = v->index;
        idxi[i] = i;
        idxj[i] = i;
    }
    
    // set indices for nonzero nondiagonal elements
    int index = (int)mesh.vertices.size();
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int i = v->index;
        
        HalfEdgeCIter he = v->he;
        do {
            int j = he->flip->vertex->index;
            
            // mosek requires only the upper (or lower triangular) part of the hessian
            if (i > j) {
                idxi[index] = i;
                idxj[index] = j;
                index++;
            }

            he = he->flip->next;
        } while (he != v->he);
    }
}

void Cetm::setEdgeLengthsAndAngles()
{
    // set edge lengths
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        lengths[e->index] = e->length() * exp(0.5*(solver.xx[he->vertex->index] +
                                                   solver.xx[he->flip->vertex->index]));
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
                if (l != -1) angles[he->index] = l == i ? M_PI : 0.0;
                i++;
                
                he = he->next;
            } while (he != f->he);
        }
    }
}

bool Cetm::computeScaleFactors()
{
    int variables = (int)mesh.vertices.size();
    int constraints = 1;
    int numanz = variables;
    
    // initialize solver
    if (!solver.initialize(variables, constraints, numanz)) return false;
    
    // setup optimization problem
    setupOptProblem();
    
    MosekSolver::MeshHandle handle((int)mesh.vertices.size(), (int)(mesh.vertices.size() + mesh.edges.size()));
    handle.computeEnergy = std::bind(&Cetm::computeEnergy, this, _1, _2);
    handle.computeGradient = std::bind(&Cetm::computeGradient, this, _1, _2);
    handle.computeHessian = std::bind(&Cetm::computeHessian, this, _1, _2);
    handle.buildGradientSparsity = std::bind(&Cetm::buildGradientSparsity, this, _1);
    handle.buildHessianSparsity = std::bind(&Cetm::buildHessianSparsity, this, _1, _2);
    solver.handle = &handle;
    
    // solve
    bool success = solver.solve(MosekSolver::GECO);
    if (success) setEdgeLengthsAndAngles();
    
    // reset solver
    solver.reset();
    
    return success;
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
