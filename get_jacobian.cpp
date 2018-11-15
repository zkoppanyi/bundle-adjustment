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
using Eigen::VectorXd;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
    // Load problem
    problem prob;    
    if (extract_problem_from_arguments(nrhs, prhs, prob) != 0) return;
        
    
    if (nlhs < 1)
    {
        LOG_ERROR("Need a return value!\n");
        return;
    }
    
    init_problem(prob);        

    // Init x0
    VectorXd x0;
    create_x0(x0, prob);
    
    // Residuals 
    VectorXd r = bundle_adjustment_fn(x0, (void*)&prob);

    // Jacobian
    Eigen::SparseMatrix<double> J(prob.sum_obs, prob.sum_unknowns); 
    bundle_adjustment_jacobian(x0, J, (void*)&prob);
       
    // Put to output
    create_sparse_matrix(&J, plhs[0]); 
    //eigen2mat(J, plhs[0]); 
    eigen2mat(r, plhs[1]); 
    
    
}
