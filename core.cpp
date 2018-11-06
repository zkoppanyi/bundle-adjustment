#include "matlab.h"

#include <vector>
#include <math.h>
#include "utils.h"
#include "core.h"
#include "problems.h"

#include <assert.h>

using namespace Eigen;
using namespace std;
        
/*
 *****************************************************
 * Calculate stochastic parameters
 *****************************************************
 */

void calc_stochastic(const optimizer_result &result, stochastic_params &params)
{
    if (result.J.rows()*result.J.cols() < MAX_MATRIX_SIZE)
    {  
        /*const MatrixXd& J = result.J;
        const VectorXd& r = result.r;

        int f = (int)J.rows() - (int)J.cols();*/

        /* 
        // calcualte rank for this; too slow...
        MatrixXd N = J.transpose()*J;
        VectorXd dr = J.transpose()*r;
        MatrixXd Naux = N;
        Naux.conservativeResize(Naux.rows(), Naux.cols()+1);
        Naux.col(Naux.cols()-1) = dr;
        FullPivLU<MatrixXd> lu_decomp(J);
        size_t rank = lu_decomp.rank();
        LOG2I("rank= ", rank);
        LOG2I("f= ", f);
        LOG2I("J.rows()= ", (int)J.rows());
        LOG2I("J.cols()= ", (int)J.cols());
        print(J);*/

        /*double sigma_0_sqr = r.dot(r) / (double)f;

        // get inverse of the normal matrix with Cheolesky
        MatrixXd N = J.transpose()*J;
        MatrixXd Ninv = N.ldlt().solve(MatrixXd::Identity(N.rows(), N.cols()));
        //MatrixXd Ninv = N.inverse();

        // Weight coefficient matrices
        MatrixXd Qxx = Ninv;
        MatrixXd Qll = J*Qxx*J.transpose();

        // Variance-covariance matrices
        params.sigma_0 = sqrt(sigma_0_sqr);
        params.Mxx = sigma_0_sqr * Ninv;
        params.Mll = sigma_0_sqr * Qll; */
        
        params.status = 0;
        strcpy(params.status_msg, "OK");
    }
    else
    {
        params.status = -1;
        sprintf(params.status_msg, "Jacobian is too large to calculate stochastic parameters!\nJacobian size: %d x %d (%d = %.1f MB)\nMax allowed: MAX_MATRIX_SIZE %d = %.1f MB", result.J.rows(), result.J.cols(), result.J.rows()*result.J.cols(), result.J.rows()*result.J.cols()*sizeof(double)/1000000.0, MAX_MATRIX_SIZE, MAX_MATRIX_SIZE*sizeof(double)/1000000.0);

    }
    
}

/*
 *****************************************************
 * Bundle adjustment
 *****************************************************
 */

int bundle_adjustment(problem &prob, problem_result &result)
{
    const double TolX = 1e-6;
	const double TolY = 1e-8;
    
    // initial guess
    VectorXd x0 = VectorXd::Zero(prob.sum_unknowns);

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

    optimizer_result opt_result;
    VectorXd sol_vec = levenberg_marquardt(bundle_adjustment_fn, bundle_adjustment_jacobian, (void*)&prob, x0, TolX, TolY, opt_result);
    //VectorXd sol_vec = simulated_annealing(bundle_adjustment_fn, bundle_adjustment_jacobian, (void*)&prob, x0, TolX, TolY, opt_result);
    result.optimizer_result = opt_result;
    
    for (size_t i = 0; i < prob.obj_pts.size(); i++)
    {
        object_pt& pt = prob.obj_pts[i];
        if (pt.type == KNOWN) continue;            
        size_t idx_pt = pt.xid;
        
        pt.x = sol_vec(idx_pt + 0);
        pt.y = sol_vec(idx_pt + 1);
        pt.z = sol_vec(idx_pt + 2);         
    }    
    
    for (size_t i = 0; i < prob.imgs.size(); i++)
    {
        img& img = prob.imgs[i];
        if (img.type == KNOWN) continue;
        
        size_t idx_img = img.xid;
        img.x = sol_vec(idx_img + 0);
        img.y = sol_vec(idx_img + 1);
        img.z = sol_vec(idx_img + 2);
        img.omega = sol_vec(idx_img + 3);
        img.phi = sol_vec(idx_img + 4);
        img.kappa = sol_vec(idx_img + 5);        
     
    }    
    
    for (size_t i = 0; i < prob.cams.size(); i++)
    {
        camera& cam = prob.cams[i];   
        if (cam.type == KNOWN) continue;
        
        size_t idx_cam = cam.xid;
        cam.f = sol_vec(idx_cam + 0);
        cam.cx = sol_vec(idx_cam + 1);
        cam.cy = sol_vec(idx_cam + 2);
        
        if (cam.cam_type == CAM_TYPE_DISTORTED)
        {
            cam.k1 = sol_vec(idx_cam + 3);
            cam.k2 = sol_vec(idx_cam + 4);
            cam.k3 = sol_vec(idx_cam + 5);
            cam.p1 = sol_vec(idx_cam + 6);
            cam.p2 = sol_vec(idx_cam + 7);           
        }
    }  
    
    for (size_t i = 0; i < prob.img_pts.size(); i++)
    {
        img img = *(prob.img_pts[i].img_ptr);
        camera cam = *( img.cam_ptr );     
        object_pt obj_pt = *( prob.img_pts[i].obj_pts_ptr );
        
        img_pt pti;
        backproject(obj_pt, img, cam, pti);
        
        prob.img_pts[i].x   = pti.x;
        prob.img_pts[i].y   = pti.y;
    }
        
        
	double r_norm = opt_result.r.norm();        
    result.residual = r_norm;    

    return 0;
}

VectorXd bundle_adjustment_fn(VectorXd x, void* params)
{
	const problem* prob = (problem*)params;
    
    VectorXd r = VectorXd::Zero(prob->sum_obs);
    
    for (size_t i = 0; i < prob->img_pts.size(); i++)
    {
        img img = *(prob->img_pts[i].img_ptr);
        camera cam = *(img.cam_ptr);     
        object_pt obj_pt = *( prob->img_pts[i].obj_pts_ptr );
        
        if (img.type == UNKNOWN) 
        {
            size_t idx_img = img.xid;
            img.x     = x(idx_img  + 0);
            img.y     = x(idx_img  + 1);
            img.z     = x(idx_img  + 2);
            img.omega = x(idx_img  + 3);
            img.phi   = x(idx_img  + 4);
            img.kappa = x(idx_img  + 5);
        }
        
        if (cam.type == UNKNOWN) 
        {        
            size_t idx_cam = cam.xid;
            cam.f     = x(idx_cam + 0);
            cam.cx    = x(idx_cam + 1);
            cam.cy    = x(idx_cam + 2);                 

            if (cam.cam_type == CAM_TYPE_DISTORTED)
            {
                cam.k1 = x(idx_cam + 3);
                cam.k2 = x(idx_cam + 4);
                cam.k3 = x(idx_cam + 5);
                cam.p1 = x(idx_cam + 6);
                cam.p2 = x(idx_cam + 7);
            }     
        }
        
        if (obj_pt.type == UNKNOWN) 
        {        
            size_t idx_pt = obj_pt.xid;
            obj_pt.x    = x(idx_pt + 0);
            obj_pt.y    = x(idx_pt + 1);
            obj_pt.z    = x(idx_pt + 2);                     
        }                
        
        img_pt pti0;
        backproject(obj_pt, img, cam, pti0);
        
        r(i*2)   = (pti0.x - prob->img_pts[i].x);
        r(i*2+1) = (pti0.y - prob->img_pts[i].y);

    }
    
    //print(r);

    return r;
}

int bundle_adjustment_jacobian(VectorXd x, Eigen::SparseMatrix<double> &J, void* params)
{
	const problem* prob = (problem*)params;
	
    //MatrixXd J = MatrixXd::Zero(prob->sum_obs, prob->sum_unknowns);
    //Eigen::SparseMatrix<double> J(prob->sum_obs, prob->sum_unknowns);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(prob->sum_obs*10);
        
    for (size_t i = 0; i < prob->img_pts.size(); i++)
    {
        img img = *(prob->img_pts[i].img_ptr);
        camera cam = *(img.cam_ptr);   
        object_pt obj_pt = *( prob->img_pts[i].obj_pts_ptr);
        
        if (img.type == UNKNOWN) 
        {
            size_t idx_img = img.xid;
            img.x     = x(idx_img  + 0);
            img.y     = x(idx_img  + 1);
            img.z     = x(idx_img  + 2);
            img.omega = x(idx_img  + 3);
            img.phi   = x(idx_img  + 4);
            img.kappa = x(idx_img  + 5);
            
            photo_jacobian_problem_exterior::populate(tripletList, (int)i*2, (int)idx_img, obj_pt, img, cam);
        }
        
        if (cam.type == UNKNOWN) 
        {
            size_t idx_cam = cam.xid;
            cam.f     = x(idx_cam + 0);
            cam.cx    = x(idx_cam + 1);
            cam.cy    = x(idx_cam + 2);

            if (cam.cam_type == CAM_TYPE_DISTORTED)
            {
                cam.k1 = x(idx_cam + 3);
                cam.k2 = x(idx_cam + 4);
                cam.k3 = x(idx_cam + 5);
                cam.p1 = x(idx_cam + 6);
                cam.p2 = x(idx_cam + 7);
                
                img_pt pt = prob->img_pts[i];
                
                //scaling for robustness
                pt.x *= 1e3;
                pt.y *= 1e3;
                photo_jacobian_problem_camera_distort::populate(tripletList, (int)i*2, (int)idx_cam, obj_pt, img, cam, pt);
            }     
            else
            {
                photo_jacobian_problem_camera::populate(tripletList, (int)i*2, (int)idx_cam, obj_pt, img, cam);
            }
        }
        
        if (obj_pt.type == UNKNOWN)  
        {        
            size_t idx_pt = obj_pt.xid;
            obj_pt.x    = x(idx_pt + 0);
            obj_pt.y    = x(idx_pt + 1);
            obj_pt.z    = x(idx_pt + 2);
            photo_jacobian_problem_triang::populate(tripletList, (int)i*2, (int)idx_pt, obj_pt, img, cam);                               
        }     
        
    }
    J.setFromTriplets(tripletList.begin(), tripletList.end());
    
    //print(J);
	return 0;
}

/*
 *****************************************************
 * Backprojection
 *****************************************************
 */

int backproject(object_pt &obj_pt, img &img, img_pt &pti)
{
    return backproject(obj_pt, img, *(img.cam_ptr), pti);
}

int backproject(const point3d &pt, img &img, img_pt &pti)
{
    return backproject(pt, img, *(img.cam_ptr), pti);
}

int backproject(object_pt &obj_pt, img &img, const camera &cam, img_pt &pti)
{
    point3d pt3d;
    pt3d.x = obj_pt.x;
    pt3d.y = obj_pt.y;
    pt3d.z = obj_pt.z;    
    pti.obj_pts_ptr = &obj_pt;    
    return backproject(pt3d, img, cam, pti);
}

int backproject(const point3d &pt, img &img, const camera &cam, img_pt &pti)
{
    double R11 = cos(img.phi)*cos(img.kappa);
    double R12 = -cos(img.phi)*sin(img.kappa);
    double R13 = sin(img.phi);
    double R21 = cos(img.omega)*sin(img.kappa)+sin(img.omega)*sin(img.phi)*cos(img.kappa);
    double R22 = cos(img.omega)*cos(img.kappa)-sin(img.omega)*sin(img.phi)*sin(img.kappa);
    double R23 = -sin(img.omega)*cos(img.phi);
    double R31 = sin(img.omega)*sin(img.kappa)-cos(img.omega)*sin(img.phi)*cos(img.kappa);
    double R32 = sin(img.omega)*cos(img.kappa)+cos(img.omega)*sin(img.phi)*sin(img.kappa);
    double R33 = cos(img.omega)*cos(img.phi);
    
    pti.x = -cam.f * (R11*(pt.x-img.x) + R12*(pt.y-img.y) + R13*(pt.z-img.z)) / (R31*(pt.x-img.x) + R32*(pt.y-img.y) + R33*(pt.z-img.z)) + cam.cx;
    pti.y = -cam.f * (R21*(pt.x-img.x) + R22*(pt.y-img.y) + R23*(pt.z-img.z)) / (R31*(pt.x-img.x) + R32*(pt.y-img.y) + R33*(pt.z-img.z)) + cam.cy;
    pti.img_ptr = &img;
    
    if (cam.cam_type == CAM_TYPE_DISTORTED)
    {   
        double x_hat = (pti.x - cam.cx) * 1e+3;
        double y_hat = (pti.y - cam.cy) * 1e+3;
        double r = sqrt(pow(x_hat, 2) + pow(y_hat, 2));
        double dx = x_hat * (cam.k1 * pow(r, 2) + cam.k2 * pow(r, 4) + cam.k3 * pow(r, 6)) + cam.p1*(pow(r,2) + 2*pow(x_hat,2)) + 2*cam.p2*x_hat*y_hat;
        double dy = y_hat * (cam.k1 * pow(r, 2) + cam.k2 * pow(r, 4) + cam.k3 * pow(r, 6)) + 2*cam.p1*x_hat*y_hat + cam.p2*(pow(r,2) + 2*pow(y_hat,2));
        
        pti.x -= dx;
        pti.y -= dy;
    }
    
    return 0;
}

int init_problem(problem &prob)
{
    size_t xid = 0;

    prob.n_unknown_imgs = 0;
    prob.start_idx_imgs = xid;
    for (size_t i = 0; i < prob.imgs.size(); i++)
    {
        if (prob.imgs[i].type == UNKNOWN) 
        {
            prob.imgs[i].xid = xid;
            prob.n_unknown_imgs++;
            xid += 6;
        }
    }   
    prob.end_idx_imgs = xid;
         
    prob.n_unknown_cams = 0;
    prob.start_idx_cams = xid;
    for (size_t i = 0; i < prob.cams.size(); i++)
    {
        if (prob.cams[i].type == UNKNOWN) 
        {
            if (prob.cams[i].cam_type == CAM_TYPE_DISTORTED)
            {
                prob.cams[i].xid = xid;
                xid += 8;
            }
            else
            {
                prob.cams[i].xid = xid;
                xid += 3;
            }   
            prob.n_unknown_cams++;
        }
    }     
    prob.end_idx_cams = xid;
    
    prob.n_unknown_obj_pts = 0;
    prob.start_idx_obj_pts = xid;
    for (size_t i = 0; i < prob.obj_pts.size(); i++)
    {
        if (prob.obj_pts[i].type == UNKNOWN) 
        {
            prob.obj_pts[i].xid = xid;
            prob.n_unknown_obj_pts++;
            xid += 3;
        }
    } 
    prob.end_idx_obj_pts = xid;
    
    prob.n_unknown_img_pts = 0;
    for (size_t i = 0; i < prob.img_pts.size(); i++)
    {
        if (prob.img_pts[i].type == UNKNOWN) prob.n_unknown_img_pts++;
    }   
            
    prob.sum_unknowns = xid;
    prob.sum_obs = prob.img_pts.size()*2;
            
    return 0;
};
















/*
 *****************************************************
 * Multi-ray triangulation
 *****************************************************
 */

/*
VectorXd multiray_triangulate_fn(VectorXd x, void* params);
MatrixXd multiray_triangulate_jacobian(VectorXd x, void* params);


int multiray_triangulate(const problem &prob, point3d &sol, problem_result &result)
{
    const double TolX = 1e-4;
	const double TolY = 1e-7;
    
    size_t n = prob.img_pts.size();
    assert(prob.obj_pts.size() > 0);
                
    // initial guess
    VectorXd x0(3);
    x0(0) = prob.obj_pts[0].x;
    x0(1) = prob.obj_pts[0].y;
    x0(2) = prob.obj_pts[0].z;
    
    optimizer_result opt_result;
    VectorXd sol_vec = levenberg_marquardt(multiray_triangulate_fn, multiray_triangulate_jacobian, (void*)&prob, x0, TolX, TolY, opt_result);
    result.optimizer_result = opt_result;
            
    sol.x = sol_vec(0);
	sol.y = sol_vec(1);
	sol.z = sol_vec(2);       

    MatrixXd J = opt_result.J;
    VectorXd r = opt_result.r;
    
	double r_norm = r.norm();        
    result.residual = r_norm;    
        
    return 0;
}

VectorXd multiray_triangulate_fn(VectorXd x, void* params)
{
	const problem* prob = (problem*)params;
	
    Eigen::VectorXd r(2*prob->img_pts.size());
    for (size_t i = 0; i < prob->img_pts.size(); i++)
    {
        img img = *(prob->img_pts[i].img_ptr);
        camera cam = *(img.cam_ptr);   
        object_pt obj_pts = *( prob->img_pts[i].obj_pts_ptr);
        
        point3d pt0;
        pt0.x = x(0);
        pt0.y = x(1);
        pt0.z = x(2);
        
        img_pt pti0;
        backproject(pt0, img, cam, pti0);
        
        r(i*2)   = (pti0.x - prob->img_pts[i].x)/4.87e-6;
        r(i*2+1) = (pti0.y - prob->img_pts[i].y)/4.87e-6;
    }
    
    //print(r);
	return r;
}

MatrixXd multiray_triangulate_jacobian(VectorXd x, void* params)
{
	const problem* prob = (problem*)params;
	
    Eigen::MatrixXd J(2*prob->img_pts.size(), 3);
    for (size_t i = 0; i < prob->img_pts.size(); i++)
    {
        img img = *(prob->img_pts[i].img_ptr);
        camera cam = *(img.cam_ptr);   
        
        photo_jacobian_problem_triang::populate(J, i*2, 0, x, img,  cam);
    }
    
    //print(J);
	return J;
}*/

