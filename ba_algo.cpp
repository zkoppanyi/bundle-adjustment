#include "matlab.h"

#include "utils.h"

#include <Eigen/Dense>
#include <iostream>
#include <string.h>

#include "core.h"
#include "matlab.h"

#include "time.h"

using namespace std;
using Eigen::MatrixXd;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
    LOG("\n\n Bundle Adjustment");
    
    // Load problem
    problem prob;    
    if (extract_problem_from_arguments(nrhs, prhs, prob) != 0) return;
        
    
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
    LOG(" ");
    LOG2I("   sum_unknowns: ", prob.sum_unknowns);
    LOG2I("        sum_obs: ", prob.sum_obs);
    LOG2I("        n_imgs : ", prob.n_imgs);
    LOG2I("  no_of_img_var: ", prob.no_of_img_var);
    LOG2I("       idx_imgs: ", prob.idx_imgs);
    LOG2I("        n_cams : ", prob.n_cams);
    LOG2I("  no_of_cam_var: ", prob.no_of_cam_var);
    LOG2I("       idx_cams: ", prob.idx_cams);

    clock_t t = clock();
    problem_result result;
    bundle_adjustment(prob, result);
    
    LOG2F("Runnning time [s]: ", (float)(clock() - t)/CLOCKS_PER_SEC);
       
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
    
    // populate returns
    create_problem_struct(prob, plhs[0]);
    
    if (nlhs > 1)
    {
        // calc stochastic parameters
        stochastic_params stoch;
        calc_stochastic(result.optimizer_result, stoch);
        
        LOG2F("Simga null: ", stoch.sigma_0);
        create_stoch_struct(result.optimizer_result, stoch, plhs[1]);
        
    }
    
}
