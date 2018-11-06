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
    VectorXd x0;
    size_t mode_idx = 2;
    
    if (extract_problem_from_arguments(nrhs, prhs, prob) != 0) return;
    init_problem(prob);        
    
    if ((nrhs == 5) || (nrhs == 4))
    {        
        /*LOG2I(" Unknown image points: ", prob.n_unknown_img_pts);
        LOG2I("       Unknown images: ", prob.n_unknown_imgs);
        LOG2I("Unknown object points: ", prob.n_unknown_obj_pts);
        LOG2I("      Unknown cameras: ", prob.n_unknown_cams);
        LOG(" ");
        LOG2I("   sum_unknowns: ", prob.sum_unknowns);
        LOG2I("        sum_obs: ", prob.sum_obs);*/

        x0 = VectorXd::Zero(prob.sum_unknowns);

        for (size_t i = 0; i < prob.obj_pts.size(); i++)
        {
            object_pt pt = prob.obj_pts[i];
            if (pt.type == KNOWN) continue;

            size_t idx_pt = pt.xid;
            x0(idx_pt + 0) = pt.x;
            x0(idx_pt + 1) = pt.y;
            x0(idx_pt + 2) = pt.z;       

        }    

        for (size_t i = 0; i < prob.imgs.size(); i++)
        {
            img img = prob.imgs[i];
            if (img.type == KNOWN) continue;

            size_t idx_img = img.xid;
            x0(idx_img + 0) = img.x;
            x0(idx_img + 1) = img.y;
            x0(idx_img + 2) = img.z;
            x0(idx_img + 3) = img.omega;
            x0(idx_img + 4) = img.phi;
            x0(idx_img + 5) = img.kappa;        

        }    

        for (size_t i = 0; i < prob.cams.size(); i++)
        {
            camera cam = prob.cams[i];           
            if (cam.type == KNOWN) continue;

            size_t idx_cam = cam.xid;
            x0(idx_cam + 0) = cam.f;
            x0(idx_cam + 1) = cam.cx;
            x0(idx_cam + 2) = cam.cy;

            if (cam.cam_type == CAM_TYPE_DISTORTED)
            {
                x0(idx_cam + 3) = cam.k1;
                x0(idx_cam + 4) = cam.k2;
                x0(idx_cam + 5) = cam.k3;
                x0(idx_cam + 6) = cam.p1;
                x0(idx_cam + 7) = cam.p2;       
            }
        }         
        
        mode_idx = 4;
        if (nrhs == 4) mode_idx = -1;
    }
    else
    {
        size_t n_x0 = mxGetM(prhs[4]);
        x0 = VectorXd::Zero(n_x0);
        double* x0_ptr = mxGetPr(prhs[4]);
        
        for (size_t i = 0; i < n_x0; i++)
        {
            x0(i) = GET(x0_ptr, i, 0, n_x0);
        }
            
        mode_idx = 5;
    }
    
    if (nlhs < 1)
    {
        LOG_ERROR("Need a return value!\n");
        return;
    }
    
    int mode = 1;
    if (nrhs >= mode_idx)
    {
        //if ( !mxIsDouble(prhs[5]) )
        {
            mode = (int)mxGetScalar(prhs[mode_idx]);
            //LOG2I("Mode: ", mode);
        }      
        //else
        //{
        //    LOG_ERROR("Mode is not a scalar value\n");
        //    return;
        //}
    }        
    
    if (mode == 1)
    {
        VectorXd r = bundle_adjustment_fn(x0, &prob);
        
        size_t n_ret = r.size();
        plhs[0] = mxCreateDoubleMatrix(n_ret, 1, mxREAL);
        double* ret_ptr = mxGetPr(plhs[0]);
        for (size_t i = 0; i < n_ret ; i++)
        {
            GET(ret_ptr, i, 0, n_ret) = r(i);            
        }           

        Eigen::SparseMatrix<double> J(prob.sum_obs, prob.sum_unknowns);  
        bundle_adjustment_jacobian(x0, J, &prob);   
        create_sparse_matrix(&J, plhs[1]); 
    }
    else if (mode == 2)
    {
        size_t n_ret = x0.size();
        plhs[0] = mxCreateDoubleMatrix(n_ret, 1, mxREAL);
        double* ret_ptr = mxGetPr(plhs[0]);
        for (size_t i = 0; i < n_ret ; i++)
        {
            GET(ret_ptr, i, 0, n_ret) = x0(i);            
        }           
    }
    else
    {
        LOG_ERROR("Mode is invalid!");
    }
    
    
    
}
