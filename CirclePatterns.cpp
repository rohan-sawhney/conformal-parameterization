#include "CirclePatterns.h"
#include "Utils.h"
#define EPSILON 1e-3

CirclePatterns::CirclePatterns(Mesh& mesh0):
Parameterization(mesh0),
angles(mesh.halfEdges.size()),
thetas(mesh.edges.size()),
radii(mesh.faces.size()-1),
eIntIndices(mesh.edges.size()),
imaginaryHe(0)
{
    
}

void CirclePatterns::setupAngleOptProblem()
{
    // set triangle sum constraint
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fIdx = f->index;
            solver.bkc[fIdx] = MSK_BK_FX;
            solver.blc[fIdx] = M_PI;
            solver.buc[fIdx] = M_PI;
        }
    }
    
    // set vertex sum constraint
    int shift = (int)(mesh.faces.size() - mesh.boundaries.size());
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        int vIdx = v->index + shift;
        if (v->isBoundary()) {
            solver.bkc[vIdx] = MSK_BK_RA;
            solver.blc[vIdx] = 0.0;
            solver.buc[vIdx] = 2*M_PI;
        
        } else {
            solver.bkc[vIdx] = MSK_BK_FX;
            solver.blc[vIdx] = 2*M_PI;
            solver.buc[vIdx] = 2*M_PI;
        }
    }
    
    // set local delaunay constraint
    shift += (int)mesh.vertices.size();
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) {
            int eIdx = eIntIndices[e->index] + shift;
            solver.bkc[eIdx] = MSK_BK_RA;
            solver.blc[eIdx] = EPSILON;
            solver.buc[eIdx] = M_PI - EPSILON;
        }
    }
    
    // set angle matrices
    int a = 0;
    int numanz = 0;
    for (HalfEdgeCIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++) {
        if (!he->onBoundary) {
            double angle = he->angle();
            
            // subtract angle from constraints and set linear constraint matrix
            solver.ptrb[a] = numanz;
            
            int fIdx = he->face->index;
            solver.blc[fIdx] -= angle;
            solver.buc[fIdx] -= angle;
            solver.asub[numanz] = fIdx;
            solver.aval[numanz] = 1.0;
            numanz++;
            
            int vIdx = he->next->next->vertex->index + (int)(mesh.faces.size() - mesh.boundaries.size());
            solver.blc[vIdx] -= angle;
            solver.buc[vIdx] -= angle;
            solver.asub[numanz] = vIdx;
            solver.aval[numanz] = 1.0;
            numanz++;
            
            if (!he->edge->isBoundary()) {
                int eIdx = eIntIndices[he->edge->index] +
                           (int)(mesh.faces.size() + mesh.vertices.size() - mesh.boundaries.size());
                solver.blc[eIdx] -= angle;
                solver.buc[eIdx] -= angle;
                solver.asub[numanz] = eIdx;
                solver.aval[numanz] = 1.0;
                numanz++;
            }
            
            solver.ptre[a] = numanz;
            
            // ensure positive and non reflex angles
            solver.bkx[a] = MSK_BK_RA;
            solver.blx[a] = EPSILON - angle;
            solver.bux[a] = (M_PI - EPSILON) - angle;
            
            // set objective
            solver.qsubi[a] = a;
            solver.qsubj[a] = a;
            solver.qval[a] = 1.0;
            solver.c[a] = 0.0;
            
            a++;
        }
    }
}

void CirclePatterns::setThetas()
{
    // set opposite angles
    int a = 0;
    for (HalfEdgeCIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++) {
        if (!he->onBoundary) angles[he->index] = he->angle() + solver.xx[a++];
        else angles[he->index] = 0.0;
    }
    
    // set thetas
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter he = e->he;
        thetas[e->index] = M_PI - angles[he->index] - angles[he->flip->index];
    }
}

bool CirclePatterns::computeAngles()
{
    int variables = 3 * (int)(mesh.faces.size() - mesh.boundaries.size()); 
    int constraints = (int)(mesh.vertices.size() + mesh.edges.size() + mesh.faces.size() -
                            imaginaryHe - mesh.boundaries.size());
    int numanz = 3 * variables - imaginaryHe;
    int numqnz = variables;
    
    // initialize solver
    if (!solver.initialize(variables, constraints, numanz, numqnz)) return false;
    
    // setup optimization problem
    setupAngleOptProblem();
    
    // solve
    bool success = solver.solve(QO);
    if (success) setThetas();
    
    // reset solver
    solver.reset();
    
    return success;
}

void CirclePatterns::setupRadiiOptProblem()
{
    // set constraint: sum of rhos over all faces equals 0
    solver.bkc[0] = MSK_BK_FX;
    solver.blc[0] = 0.0;
    solver.buc[0] = 0.0;
    
    // set linear constraint vector and bounds
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fIdx = f->index;
            
            solver.ptrb[fIdx] = fIdx;
            solver.ptre[fIdx] = fIdx+1;
            solver.asub[fIdx] = 0;
            solver.aval[fIdx] = 1.0;
            
            solver.bkx[fIdx] = MSK_BK_FR;
            solver.blx[fIdx] = -MSK_INFINITY;
            solver.bux[fIdx] = MSK_INFINITY;
            
            solver.c[fIdx] = 1.0;
        }
    }
}

double ImLi2Sum(double dp, double theta)
{
    double tStar = M_PI - theta;
    double x = 2*atan(tanh(0.5*dp) * tan(0.5*tStar));
    
    return x*dp + Cl2(x + tStar) + Cl2(-x + tStar) - Cl2(2.0*tStar);
}

double fe(double dp, double theta)
{
    return atan2(sin(theta), exp(dp) - cos(theta));
}

void CirclePatterns::computeEnergy(double& energy, const double *rho)
{
    energy = 0.0;
    
    // sum over edges
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        int fk = e->he->face->index;
        
        if (e->isBoundary()) {
            energy -= 2*(M_PI - thetas[e->index]) * rho[fk];
            
        } else {
            int fl = e->he->flip->face->index;
            energy += ImLi2Sum(rho[fk] - rho[fl], thetas[e->index]) -
                      (M_PI - thetas[e->index])*(rho[fk] + rho[fl]);
        }
    }
    
    // sum over faces
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) energy += 2*M_PI*rho[f->index];
    }
}

void CirclePatterns::computeGradient(double *gradient, const double *rho)
{
    // loop over faces
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fk = f->index;
            gradient[fk] = 2*M_PI;
            
            // sum of adjacent edges
            HalfEdgeCIter he = f->he;
            do {
                EdgeCIter e = he->edge;
                if (e->isBoundary()) {
                    gradient[fk] -= 2*(M_PI - thetas[e->index]);
                    
                } else {
                    HalfEdgeCIter h = e->he;
                    int fl = fk == (int)h->face->index ? h->flip->face->index : h->face->index;
                    gradient[fk] -= 2*fe(rho[fk] - rho[fl], thetas[e->index]);
                }
                
                he = he->next;
            } while (he != f->he);
        }
    }
}

void CirclePatterns::computeHessian(double *hessian, const double *rho)
{
    // sum over faces
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) hessian[f->index] = 0.0;
    }
    
    // sum over edges
    int shift = (int)(mesh.faces.size() - mesh.boundaries.size());
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) {
            int fk = e->he->face->index;
            int fl = e->he->flip->face->index;
            
            // mosek requires only the upper (or lower triangular) part of the hessian
            if (fk < fl) std::swap(fk, fl);
            
            double hessval = sin(thetas[e->index]) / (cosh(rho[fk] - rho[fl]) - cos(thetas[e->index]));
            hessian[fk] += hessval;
            hessian[fl] += hessval; // cosh is symmetric
            hessian[shift + eIntIndices[e->index]] = -hessval;
        }
    }
}

void CirclePatterns::buildGradientSparsity(int *idx)
{
    // set indices for nonzero elements
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fIdx = f->index;
            idx[fIdx] = fIdx;
        }
    }
}

void CirclePatterns::buildHessianSparsity(int *idxi, int *idxj)
{
    // set indices for nonzero diagonal elements
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int fIdx = f->index;
            idxi[fIdx] = fIdx;
            idxj[fIdx] = fIdx;
        }
    }
    
    // set indices for nonzero nondiagonal elements
    int shift = (int)(mesh.faces.size() - mesh.boundaries.size());
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) {
            int fk = e->he->face->index;
            int fl = e->he->flip->face->index;
            
            // mosek requires only the upper (or lower triangular) part of the hessian
            if (fk < fl) std::swap(fk, fl);
            
            idxi[shift + eIntIndices[e->index]] = fk;
            idxj[shift + eIntIndices[e->index]] = fl;
        }
    }
}

void CirclePatterns::setRadii()
{
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) radii[f->index] = exp(solver.xx[f->index]);
    }
}

bool CirclePatterns::computeRadii()
{
    int variables = (int)(mesh.faces.size() - mesh.boundaries.size());
    int constraints = 1;
    int numanz = variables;
    
    // initialize solver
    if (!solver.initialize(variables, constraints, numanz)) return false;
    
    // setup optimization problem
    setupRadiiOptProblem();
    
    MeshHandle handle((int)(mesh.faces.size()-mesh.boundaries.size()),
                      (int)(mesh.faces.size()+mesh.edges.size()-imaginaryHe-mesh.boundaries.size()));
    handle.computeEnergy = std::bind(&CirclePatterns::computeEnergy, this, _1, _2);
    handle.computeGradient = std::bind(&CirclePatterns::computeGradient, this, _1, _2);
    handle.computeHessian = std::bind(&CirclePatterns::computeHessian, this, _1, _2);
    handle.buildGradientSparsity = std::bind(&CirclePatterns::buildGradientSparsity, this, _1);
    handle.buildHessianSparsity = std::bind(&CirclePatterns::buildHessianSparsity, this, _1, _2);
    solver.handle = &handle;
    
    // solve
    bool success = solver.solve(GECO);
    if (success) setRadii();
    
    // reset solver
    solver.reset();
    
    return success;
}

void CirclePatterns::computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths)
{
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        HalfEdgeCIter h1 = e->he;
        
        if (e->isBoundary()) {
            angles[h1->index] = M_PI - thetas[e->index];
        
        } else {
            HalfEdgeCIter h2 = h1->flip;
            double dp = log(radii[h1->face->index]) - log(radii[h2->face->index]);
            angles[h1->index] = fe(dp, thetas[e->index]);
            angles[h2->index] = fe(-dp, thetas[e->index]);
        }
        
        lengths[e->index] = 2.0*radii[h1->face->index]*sin(angles[h1->index]);
    }
}

void CirclePatterns::performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir,
                                       Eigen::VectorXd& lengths, std::unordered_map<int, bool>& visited,
                                       std::stack<EdgeCIter>& stack)
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

void CirclePatterns::setUVs()
{
    // compute edge lengths
    Eigen::VectorXd lengths(mesh.edges.size());
    computeAnglesAndEdgeLengths(lengths);
    
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
        
        performFaceLayout(h1, dir, lengths, visited, stack);
        performFaceLayout(h2, -dir, lengths, visited, stack);
    }
    
    normalize();
}

void CirclePatterns::parameterize()
{
    // set interior edge indices
    int eIdx = 0;
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!e->isBoundary()) eIntIndices[e->index] = eIdx++;
        else {
            eIntIndices[e->index] = -1;
            imaginaryHe++;
        }
    }
    
    // compute angles
    if (!computeAngles()) {
        std::cout << "Unable to compute angles" << std::endl;
        return;
    }
    
    // compute radii
    if (!computeRadii()) {
        std::cout << "Unable to compute radii" << std::endl;
        return;
    }
    
    // set uvs
    setUVs();
}
