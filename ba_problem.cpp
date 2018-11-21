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
                
    init_problem(prob);        

    // Init x0
    VectorXd x0;
    if (nrhs > 4)
    {
        mat2eigen(prhs[4], x0 );
    }
    else
    {        
        create_x0(x0, prob);
    }
    
    // Residuals 
    VectorXd r = bundle_adjustment_fn(x0, (void*)&prob);

    // Jacobian
    Eigen::SparseMatrix<double> J(prob.sum_obs, prob.sum_unknowns); 
    bundle_adjustment_jacobian(x0, J, (void*)&prob);
       
    // Put to output
    create_sparse_matrix(&J, plhs[0]); 
    //eigen2mat(J, plhs[0]); 
    eigen2mat(r, plhs[1]); 
    
    //get indeces
    if (nlhs > 2)
    {
        size_t n_ret = prob.obj_pts.size();
        mxArray* m_idx_obj = mxCreateDoubleMatrix(n_ret, 3, mxREAL);
        double* idx_obj = mxGetPr(m_idx_obj);    
        for (size_t i = 0; i < prob.obj_pts.size(); i++)
        {
            object_pt pt = prob.obj_pts[i];
            if (pt.type == KNOWN)
            {
                GET(idx_obj, i, 0, n_ret) = pt.id;
                GET(idx_obj, i, 1, n_ret) = -1;  
                GET(idx_obj, i, 2, n_ret) = -1;
            }
            else
            {
                GET(idx_obj, i, 0, n_ret) = pt.id;
                GET(idx_obj, i, 1, n_ret) = pt.xid+1;
                GET(idx_obj, i, 2, n_ret) = pt.xid+3;
            }   
        }    

        n_ret = prob.imgs.size();
        mxArray* m_idx_imgs = mxCreateDoubleMatrix(n_ret, 3, mxREAL);
        double* idx_imgs = mxGetPr(m_idx_imgs);    
        for (size_t i = 0; i < prob.imgs.size(); i++)
        {
            img img = prob.imgs[i];
            if (img.type == KNOWN)
            {
                GET(idx_imgs, i, 0, n_ret) = img.id;
                GET(idx_imgs, i, 1, n_ret) = -1;       
                GET(idx_imgs, i, 2, n_ret) = -1;
            }
            else
            {
                GET(idx_imgs, i, 0, n_ret) = img.id;
                GET(idx_imgs, i, 1, n_ret) = img.xid+1;
                GET(idx_imgs, i, 2, n_ret) = img.xid+6;            
            }          

        }   

        n_ret = prob.cams.size();
        mxArray* m_idx_cams = mxCreateDoubleMatrix(n_ret, 3, mxREAL);
        double* idx_cams = mxGetPr(m_idx_cams);  
        for (size_t i = 0; i < prob.cams.size(); i++)
        {
            camera cam = prob.cams[i];           
            
            if (cam.type == KNOWN)
            {
                GET(idx_cams, i, 0, n_ret) = cam.id;
                GET(idx_cams, i, 1, n_ret) = -1;       
                GET(idx_cams, i, 2, n_ret) = -1;
            }
            else
            {
                GET(idx_cams, i, 0, n_ret) = cam.id;
                GET(idx_cams, i, 1, n_ret) = cam.xid+1;       
                GET(idx_cams, i, 2, n_ret) = cam.xid+3;

                if (cam.cam_type == CAM_TYPE_DISTORTED)
                {
                   GET(idx_cams, i, 2, n_ret) = cam.xid+8;       
                }
            }              
        }  
        mxArray* m_x0;
        eigen2mat(x0, m_x0);
        
        const char *fieldnames[] = {"x0", "idx_imgs", "idx_cams", "idx_obj_pts"};        
        plhs[2]  = mxCreateStructMatrix(1,1,4,fieldnames);
        mxSetFieldByNumber(plhs[2], 0,0, m_x0);
        mxSetFieldByNumber(plhs[2], 0,1, m_idx_imgs);
        mxSetFieldByNumber(plhs[2], 0,2, m_idx_cams);
        mxSetFieldByNumber(plhs[2], 0,3, m_idx_obj);
    }
    
}
