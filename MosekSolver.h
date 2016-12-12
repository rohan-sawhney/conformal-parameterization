#ifndef MOSEK_SOLVER_H
#define MOSEK_SOLVER_H

#include <mosek/h/mosek.h>
#include <functional>

namespace MosekSolver {
    struct MeshHandle {
        // typedefs
        typedef std::function<void(double&, const double*)> ComputeEnergy;
        typedef std::function<void(double*, const double*)> ComputeGradient;
        typedef std::function<void(double*, const double*)> ComputeHessian;
        typedef std::function<void(int*)> BuildGradientSparsity;
        typedef std::function<void(int*, int*)> BuildHessianSparsity;
        
        // Constructor
        MeshHandle(int numGrad_, int numHess_):
        numGrad(numGrad_),
        numHess(numHess_)
        {}
        
        // Member variables
        int numGrad;
        int numHess;
        ComputeEnergy computeEnergy;
        ComputeGradient computeGradient;
        ComputeHessian computeHessian;
        BuildGradientSparsity buildGradientSparsity;
        BuildHessianSparsity buildHessianSparsity;
    };
    
    enum ProblemType {
        QO, // quadratic
        GECO // general
    };
    
    class Solver {
    public:
        // Constructor
        Solver();
        
        // Destructor
        ~Solver();
        
        // Initialize
        bool initialize(MSKint32t variables_, MSKint32t constraints_,
                        MSKint32t numanz_, MSKint32t numqcnz_, MSKint32t numqnz_);
        
        // Solve
        bool solve(ProblemType type);
        
        // Reset
        void reset();
        
        // Member variables
        MSKboundkeye *bkc, *bkx;
        MSKint32t *ptrb, *ptre, *asub, *qcsubi, *qcsubj, *qsubi, *qsubj;
        MSKrealt *blc, *buc, *blx, *bux, *aval, *qcval, *qval, *c, *xx;
        MeshHandle *handle;
        
    private:
        // Initializes environment
        void initializeEnvironment();
        
        // Creates task
        bool createTask(MSKtask_t& task);
        
        // Allocates internal memory
        void allocateMemory();
        
        // Sets solution
        bool setSolution();
        
        // Deallocates internal memory
        void deallocateMemory();
        
        // Member variables
        MSKenv_t env;
        MSKtask_t task;
        MSKrescodee r;
        MSKint32t variables, constraints, numanz, numqcnz, numqnz;
    };
}

#endif
