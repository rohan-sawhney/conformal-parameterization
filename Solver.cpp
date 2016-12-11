#include "Solver.h"
#include <deque>
#include <eigen/SparseCholesky>
#define beta 0.9
#define EPSILON 1e-9

Solver::Solver(int n0):
n(n0)
{

}

void Solver::gradientDescent()
{
    int k = 1;
    double f = 0.0;
    x = Eigen::VectorXd::Zero(n);
    handle->computeEnergy(f, x);
    Eigen::VectorXd xp = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);

    while (true) {
        // compute momentum term
        v = x;
        if (k > 1) v += (k-2)*(x - xp)/(k+1);
        
        // compute update direction
        Eigen::VectorXd g(n);
        handle->computeGradient(g, v);
        
        // compute step size
        double t = 1.0;
        double fp = 0.0;
        Eigen::VectorXd xn = v - t*g;
        Eigen::VectorXd xnv = xn - v;
        handle->computeEnergy(fp, v);
        handle->computeEnergy(f, xn);
        while (f > fp + g.dot(xnv) + xnv.dot(xnv)/(2*t)) {
            t = beta*t;
            xn = v - t*g;
            xnv = xn - v;
            handle->computeEnergy(f, xn);
        }

        // update
        xp = x;
        x = xn;
        k++;
        
        // check termination condition
        if (fabs(f - fp) < EPSILON) break;
    }

    std::cout << "f: " << f << " k: " << k << std::endl;
}

void solve(Eigen::VectorXd& x, const Eigen::VectorXd& b, const Eigen::SparseMatrix<double>& A)
{
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A);
    x = solver.solve(b);
}

void Solver::newton()
{
    int k = 0;
    double f = 0.0;
    x = Eigen::VectorXd::Zero(n);
    handle->computeEnergy(f, x);

    const double alpha = 0.5;
    while (true) {
        // compute update direction
        Eigen::VectorXd g(n);
        handle->computeGradient(g, x);

        Eigen::SparseMatrix<double> H(n, n);
        handle->computeHessian(H, x);

        Eigen::VectorXd p;
        solve(p, g, H);

        // compute step size
        double t = 1.0;
        double fp = f;
        handle->computeEnergy(f, x - t*p);
        while (f > fp - alpha*t*g.dot(p)) {
            t = beta*t;
            handle->computeEnergy(f, x - t*p);
        }

        // update
        x -= t*p;
        k++;

        // check termination condition
        if (fabs(f - fp) < EPSILON) break;
    }
    
    std::cout << "f: " << f << " k: " << k << std::endl;
}

void Solver::lbfgs(int m)
{
    int k = 0;
    double f = 0.0;
    x = Eigen::VectorXd::Zero(n);
    handle->computeEnergy(f, x);
    Eigen::VectorXd g(n);
    handle->computeGradient(g, x);
    std::deque<Eigen::VectorXd> s;
    std::deque<Eigen::VectorXd> y;
    
    const double alpha = 1e-4;
    while (true) {
        // compute update direction
        const int l = std::min(k, m);
        Eigen::VectorXd q = -g;
        
        Eigen::VectorXd a(l);
        for (int i = l-1; i >= 0; i--) {
            a(i) = s[i].dot(q) / y[i].dot(s[i]);
            q -= a(i)*y[i];
        }
        
        Eigen::VectorXd p = q;
        if (l > 0) p *= y[l-1].dot(s[l-1]) / y[l-1].dot(y[l-1]);
        
        for (int i = 0; i < l; i++) {
            double b = y[i].dot(p) / y[i].dot(s[i]);
            p += (a(i) - b)*s[i];
        }
        
        // compute step size
        double t = 1.0;
        double fp = f;
        handle->computeEnergy(f, x + t*p);
        while (f > fp + alpha*t*g.dot(p)) {
            t = beta*t;
            handle->computeEnergy(f, x + t*p);
        }
        
        // update
        Eigen::VectorXd xp = x;
        Eigen::VectorXd gp = g;
        x += t*p;
        handle->computeGradient(g, x);
        k++;
        
        // update history
        if (k > m) {
            s.pop_front();
            y.pop_front();
        }
        s.push_back(x - xp);
        y.push_back(g - gp);
        
        // check termination condition
        if (fabs(f - fp) < EPSILON) break;
    }
    
    std::cout << "f: " << f << " k: " << k << std::endl;
}
