#include <vector>

#include "utils.h"
#include "structs.h"
#include "core.h"
#include "optim.h"

#include "matlab.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Load problem
    problem prob;    
    extract_problem_from_arguments(nrhs, prhs, prob);
    
    if (nlhs < 1)
    {
        LOG_ERROR("Need a return value!\n");
        return;
    }
    
    init_problem(prob);    
    LOG2I(" Unknown image points: ", prob.n_unknown_img_pts);
    LOG2I("       Unknown images: ", prob.n_unknown_imgs);
    LOG2I("Unknown object points: ", prob.n_unknown_obj_pts);
    LOG2I("      Unknown cameras: ", prob.n_unknown_cams);
    
    point3d pt_sol;
    problem_result result;
    multiray_triangulate(prob, pt_sol, result);
    
    if (result.optimizer_result.stopping_criteria == THRESHOLD_REACHED)
    {
        LOG("Optimizer: Threshold reached!");
    }
    else if (result.optimizer_result.stopping_criteria == MAXIMUM_ITERATION_REACHED)
    {
        LOG_WARNING("Optimizer: Maximum iteration number reached! Check results!");
    }
    else if (result.optimizer_result.stopping_criteria == UNBALANCED_PROBLEM)
    {
        LOG_WARNING("Unbalanced problem! Check results!");
    }
    else if (result.optimizer_result.stopping_criteria == ERROR)
    {
        LOG_ERROR("Optimizer reported error!");
    }   
    
    LOG2I("Optimizer: Number of iterations: ", result.optimizer_result.no_of_iterations);
    LOG2F("L2 norm of residuals: ", result.residual);
    LOG2F("L2 norm of residuals per no. of unknowns: ", result.residual / result.optimizer_result.r.size());
    
    //calc_stochastic(result.optimizer_result);
            
    plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);
    double* ret = mxGetPr(plhs[0]);
    ret[0] = pt_sol.x;
    ret[1] = pt_sol.y;
    ret[2] = pt_sol.z;

    
}