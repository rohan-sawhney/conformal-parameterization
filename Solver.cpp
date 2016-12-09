#include "Solver.h"
#include <eigen/SparseCholesky>
#define EPSILON 1e-6
#define alpha 0.2
#define beta 0.9

Solver::Solver(int n0):
n(n0),
x(Eigen::VectorXd::Zero(n0))
{

}

void Solver::gradientDescent()
{
    // TODO: acceleration
    double f = 0.0;
    handle->computeObjective(f, x);

    while (true) {
        // compute update direction
        Eigen::VectorXd g(n);
        handle->computeGradient(g, x);

        // compute step size
        double fp = f;
        double t = 1.0;
        handle->computeObjective(f, x - t*g);
        while (f > fp - alpha*t*g.squaredNorm()) {
            t = beta*t;
            handle->computeObjective(f, x - t*g);
        }

        // update
        x -= t*g;

        // check termination condition
        if (fabs(f - fp) < EPSILON) break;
    };
}

void Solver::coordinateDescent()
{
    // TODO
}

void solve(Eigen::VectorXd& x, const Eigen::VectorXd& g, const Eigen::SparseMatrix<double>& H)
{
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(H);
    x = solver.solve(g);
}

void Solver::newton()
{   
    double f = 0.0;
    handle->computeObjective(f, x);

    while (true) {
        // compute update direction
        Eigen::VectorXd g(n);
        handle->computeGradient(g, x);

        Eigen::SparseMatrix<double> H(n, n);
        handle->computeHessian(H, x);

        Eigen::VectorXd p;
        solve(p, g, H);

        // compute step size
        double fp = f;
        double t = 1.0;
        handle->computeObjective(f, x - t*p);
        while (f > fp - alpha*t*g.dot(p)) {
            t = beta*t;
            handle->computeObjective(f, x - t*p);
        }

        // update
        x -= t*p;

        // check termination condition
        if (fabs(f - fp) < EPSILON) break;
    };
}

void Solver::lbfgs()
{
    // TODO
}
