#include "MosekSolver.h"
using namespace MosekSolver;

Solver::Solver():
bkc(NULL), bkx(NULL),
ptrb(NULL), ptre(NULL), asub(NULL), qsubi(NULL), qsubj(NULL),
blc(NULL), buc(NULL), blx(NULL), bux(NULL), aval(NULL), qval(NULL), c(NULL), xx(NULL),
handle(NULL), 
env(NULL),
task(NULL),
r(MSK_RES_OK),
variables(0),
constraints(0),
numanz(0),
numqcnz(0),
numqnz(0)
{
    initializeEnvironment();
}

Solver::~Solver()
{
    MSK_deleteenv(&env);
}

static void MSKAPI printError(void *handle, MSKCONST char err[])
{
    printf("%s", err);
}

void Solver::initializeEnvironment()
{
    r = MSK_makeenv(&env, NULL);
    if (r == MSK_RES_OK) {
        MSK_linkfunctoenvstream(env, MSK_STREAM_LOG, NULL, printError);
        r = MSK_initenv(env);
    }
}

bool Solver::createTask(MSKtask_t& task)
{
    r = MSK_maketask(env, constraints, variables, &task);
    if (r != MSK_RES_OK) return false;
    
    //MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printError);
    
    return true;
}

void Solver::allocateMemory()
{
    bkc = new MSKboundkeye[constraints]; blc = new MSKrealt[constraints]; buc = new MSKrealt[constraints];
    bkx = new MSKboundkeye[variables]; blx = new MSKrealt[variables]; bux = new MSKrealt[variables];
    ptrb = new MSKint32t[variables]; ptre = new MSKint32t[variables];
    asub = new MSKint32t[numanz]; aval = new MSKrealt[numanz];
    qcsubi = new MSKint32t[numqcnz]; qcsubj = new MSKint32t[numqcnz]; qcval = new MSKrealt[numqcnz];
    qsubi = new MSKint32t[numqnz]; qsubj = new MSKint32t[numqnz];
    qval = new MSKrealt[numqnz]; c = new MSKrealt[variables],
    xx = new MSKrealt[variables];
}

bool Solver::initialize(MSKint32t variables_, MSKint32t constraints_,
                        MSKint32t numanz_, MSKint32t numqcnz_, MSKint32t numqnz_)
{
    if (r != MSK_RES_OK) return false;
    
    variables = variables_;
    constraints = constraints_;
    numanz = numanz_;
    numqcnz = numqcnz_;
    numqnz = numqnz_;
    
    if (!createTask(task)) return false;
    allocateMemory();
    
    return true;
}

static MSKint32t MSKAPI nlgetva(MSKuserhandle_t nlhandle, MSKCONST MSKrealt *xx,
                                MSKrealt yo, MSKCONST MSKrealt *yc, MSKrealt *objval,
                                MSKint32t *numgrdobjnz, MSKint32t *grdobjsub, MSKrealt *grdobjval,
                                MSKint32t numi, MSKCONST MSKint32t *subi, MSKrealt *conval,
                                MSKCONST MSKint32t *grdconptrb, MSKCONST MSKint32t *grdconptre,
                                MSKCONST MSKint32t *grdconsub, MSKrealt *grdconval, MSKrealt *grdlag,
                                MSKint32t maxnumhesnz, MSKint32t *numhesnz,
                                MSKint32t *hessubi, MSKint32t *hessubj, MSKrealt *hesval)
{
    MeshHandle *handle = (MeshHandle *)nlhandle;
    
    // Compute energy
    if (objval) handle->computeEnergy(objval[0], xx);
    
    // Set gradient sparsity
    if (numgrdobjnz) numgrdobjnz[0] = handle->numGrad;
    if (grdobjsub) handle->buildGradientSparsity(grdobjsub);
    
    if (grdobjval) handle->computeGradient(grdobjval, xx);
    
    // Set constraint val to 0
    if (conval) {
        for (MSKint32t i = 0; i < numi; i++) conval[i] = 0.0;
    }
    
    if (grdlag) grdlag = grdobjval;
    
    // Compute hessian
    if (maxnumhesnz) {
        if (yo != 0.0) {
            numhesnz[0] = handle->numHess;
            if (hessubi && hessubj && hesval) {
                handle->buildHessianSparsity(hessubi, hessubj);
                handle->computeHessian(hesval, xx);
            }
            
        } else {
            numhesnz[0] = 0;
        }
    }
    
    return MSK_RES_OK;
}

static MSKint32t MSKAPI nlgetsp(MSKuserhandle_t nlhandle, MSKint32t *numgrdobjnz, MSKint32t *grdobjsub,
                                MSKint32t i, MSKbooleant *convali, MSKint32t *grdconinz, MSKint32t *grdconisub,
                                MSKint32t yo, MSKint32t numycnz, MSKCONST MSKint32t *ycsub,
                                MSKint32t maxnumhesnz, MSKint32t *numhesnz,
                                MSKint32t *hessubi, MSKint32t *hessubj)
{
    MeshHandle *handle = (MeshHandle *)nlhandle;
    
    // Set gradient sparsity
    if (numgrdobjnz) numgrdobjnz[0] = handle->numGrad;
    if (grdobjsub) handle->buildGradientSparsity(grdobjsub);
    
    if (convali) convali[0] = 0;
    if (grdconinz) grdconinz[0] = 0;
    
    // Set hessian sparsity
    if (maxnumhesnz) {
        if (yo != 0.0) {
            numhesnz[0] = handle->numHess;
            if (hessubi && hessubj) handle->buildHessianSparsity(hessubi, hessubj);
            
        } else {
            numhesnz[0] = 0;
        }
    }
    
    return MSK_RES_OK;
}

bool Solver::setSolution()
{
    MSKsolstae solsta;
    MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
    
    switch (solsta) {
        case MSK_SOL_STA_OPTIMAL:
        case MSK_SOL_STA_NEAR_OPTIMAL:
            r = MSK_getxx(task, MSK_SOL_ITR, xx);
            break;
        default:
            return false;
    }
    
    return r == MSK_RES_OK;
}

bool Solver::solve(ProblemType type)
{
    // Set input data
    r = MSK_inputdata(task, constraints, variables, constraints, variables, c, 0,
                      ptrb, ptre, asub, aval, bkc, blc, buc, bkx, blx, bux);
    if (r != MSK_RES_OK) return false;
    
    switch (type) {
        case QO:
            r = MSK_putqobj(task, numqnz, qsubi, qsubj, qval);
            if (numqcnz > 0) r = MSK_putqconk(task, 0, numqcnz, qcsubi, qcsubj, qcval);
            break;

        case GECO:
            r = MSK_putnlfunc(task, handle, nlgetsp, nlgetva);
            break;
            
        default:
            break;
    }
    if (r != MSK_RES_OK) return false;
   
    // Optimize
    r = MSK_optimize(task);
    if (r != MSK_RES_OK) return false;
    
    // Set solution
    return setSolution();
}

void Solver::deallocateMemory()
{
    delete [] bkc; delete [] blc; delete [] buc;
    delete [] bkx; delete [] blx; delete [] bux;
    delete [] ptrb; delete [] ptre; delete [] asub; delete [] aval;
    delete [] qcsubi; delete [] qcsubj; delete [] qcval;
    delete [] qsubi; delete [] qsubj; delete [] qval; delete [] c;
    delete [] xx;
}

void Solver::reset()
{
    deallocateMemory();
    MSK_deletetask(&task);
}
