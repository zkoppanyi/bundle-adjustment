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
     const MatrixXd& J = result.J;
     const VectorXd& r = result.r;
         
     int f = (int)J.rows() - (int)J.cols();
    
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

    double sigma_0_sqr = r.dot(r) / (double)f;

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
    params.Mll = sigma_0_sqr * Qll;  
}

/*
 *****************************************************
 * Bundle adjustment
 *****************************************************
 */

VectorXd bundle_adjustment_fn(VectorXd x, void* params);
MatrixXd bundle_adjustment_jacobian(VectorXd x, void* params);

int bundle_adjustment(const problem &prob, problem &sol, problem_result &result)
{
    const double TolX = 1e-6;
	const double TolY = 1e-8;
    
    sol = prob;
    init_problem(sol);
    
    // initial guess
    VectorXd x0 = VectorXd::Zero(sol.sum_unknowns);
    for (size_t i = 0; i < sol.n_imgs; i++)
    {
        img img = sol.imgs[i];
        size_t img_id = img.id;
    
        int cam_id = img.cam_id;
        camera cam = sol.cams[cam_id]; 
        
        size_t idx_img = sol.idx_imgs + sol.no_of_img_var*img_id;
        x0(idx_img + 0) = img.x;
        x0(idx_img + 1) = img.y;
        x0(idx_img + 2) = img.z;
        x0(idx_img + 3) = img.omega;
        x0(idx_img + 4) = img.phi;
        x0(idx_img + 5) = img.kappa;
        
        size_t idx_cam = sol.idx_cams + sol.no_of_cam_var*cam_id;
        x0(idx_cam + 0) = cam.f;
        x0(idx_cam + 1) = cam.cx;
        x0(idx_cam + 2) = cam.cy;
        
        if (cam.cam_type == CAM_TYPE_DISTORTED)
        {
            x0(idx_cam + 3) = cam.k1 / 1e-4;
            x0(idx_cam + 4) = cam.k2 / 1e-6;
            x0(idx_cam + 5) = cam.k3 / 1e-10;
            x0(idx_cam + 6) = cam.p1 / 1e-5;
            x0(idx_cam + 7) = cam.p2 / 1e-4;           
        }
    }    
    
    optimizer_result opt_result;
    VectorXd sol_vec = levenberg_marquardt(bundle_adjustment_fn, bundle_adjustment_jacobian, (void*)&sol, x0, TolX, TolY, opt_result);
    //VectorXd sol_vec = simulated_annealing(bundle_adjustment_fn, bundle_adjustment_jacobian, (void*)&sol, x0, TolX, TolY, opt_result);
    result.optimizer_result = opt_result;
    
    for (size_t i = 0; i < sol.n_imgs; i++) 
    {
        img &img = sol.imgs[i];
        size_t img_id = img.id;

        int cam_id = img.cam_id;
        camera& cam = sol.cams[cam_id]; 
        
        size_t idx_img = sol.idx_imgs + sol.no_of_img_var*img_id;
        img.x     = sol_vec(idx_img  + 0);
        img.y     = sol_vec(idx_img  + 1);
        img.z     = sol_vec(idx_img  + 2);
        img.omega = sol_vec(idx_img  + 3);
        img.phi   = sol_vec(idx_img  + 4);
        img.kappa = sol_vec(idx_img  + 5);
                    
        size_t idx_cam = sol.idx_cams + sol.no_of_cam_var*cam_id;
        cam.f     = sol_vec(idx_cam + 0);
        cam.cx    = sol_vec(idx_cam + 1);
        cam.cy    = sol_vec(idx_cam + 2);      
        
        if (cam.cam_type == CAM_TYPE_DISTORTED)
        {
            cam.k1 = sol_vec(idx_cam + 3) * 1e-4;
            cam.k2 = sol_vec(idx_cam + 4) * 1e-6;
            cam.k3 = sol_vec(idx_cam + 5) * 1e-10;
            cam.p1 = sol_vec(idx_cam + 6) * 1e-5;
            cam.p2 = sol_vec(idx_cam + 7) * 1e-4;
        }

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
        int img_id = prob->img_pts[i].img_id;
        img img = prob->imgs[img_id];
        
        int cam_id = img.cam_id;
        camera cam = prob->cams[cam_id];        
        
        size_t idx_img = prob->idx_imgs + prob->no_of_img_var*img_id;
        img.x     = x(idx_img  + 0);
        img.y     = x(idx_img  + 1);
        img.z     = x(idx_img  + 2);
        img.omega = x(idx_img  + 3);
        img.phi   = x(idx_img  + 4);
        img.kappa = x(idx_img  + 5);
        
        size_t idx_cam = prob->idx_cams + prob->no_of_cam_var*cam_id;
        cam.f     = x(idx_cam + 0);
        cam.cx    = x(idx_cam + 1);
        cam.cy    = x(idx_cam + 2);        
        
        if (cam.cam_type == CAM_TYPE_DISTORTED)
        {
            cam.k1 = x(idx_cam + 3) * 1e-4;
            cam.k2 = x(idx_cam + 4) * 1e-6;
            /*cam.k3 = x(idx_cam + 5) * 1e-10;
            cam.p1 = x(idx_cam + 6) * 1e-5;
            cam.p2 = x(idx_cam + 7) * 1e-4;*/
        }        
        
        int obj_pts_id = prob->img_pts[i].obj_pts_id;
        object_pt obj_pts = prob->obj_pts[obj_pts_id];
        
        img_pt pti0;
        backproject(obj_pts, img, cam, pti0);
        
        r(i*2)   = pti0.x - prob->img_pts[i].x;
        r(i*2+1) = pti0.y - prob->img_pts[i].y;
    }
    
    //print(r);
	return r;
}

MatrixXd bundle_adjustment_jacobian(VectorXd x, void* params)
{
	const problem* prob = (problem*)params;
	
    MatrixXd J = MatrixXd::Zero(prob->sum_obs, prob->sum_unknowns);
    
    for (size_t i = 0; i < prob->img_pts.size(); i++)
    {
        int img_id = prob->img_pts[i].img_id;
        img img = prob->imgs[img_id];
        
        int cam_id = img.cam_id;
        camera cam = prob->cams[cam_id];  
        
        int obj_pts_id = prob->img_pts[i].obj_pts_id;
        object_pt obj_pts = prob->obj_pts[obj_pts_id];
        
        size_t idx_img = prob->idx_imgs + prob->no_of_img_var*img_id;
        img.x     = x(idx_img  + 0);
        img.y     = x(idx_img  + 1);
        img.z     = x(idx_img  + 2);
        img.omega = x(idx_img  + 3);
        img.phi   = x(idx_img  + 4);
        img.kappa = x(idx_img  + 5);
        
        size_t idx_cam = prob->idx_cams + prob->no_of_cam_var*cam_id;
        cam.f     = x(idx_cam + 0);
        cam.cx    = x(idx_cam + 1);
        cam.cy    = x(idx_cam + 2);
        
        if (cam.cam_type == CAM_TYPE_DISTORTED)
        {
            cam.k1 = x(idx_cam + 3) * 1e-4;
            cam.k2 = x(idx_cam + 4) * 1e-6;
            /*cam.k3 = x(idx_cam + 5) * 1e-10;
            cam.p1 = x(idx_cam + 6) * 1e-5;
            cam.p2 = x(idx_cam + 7) * 1e-4;*/
        }        
        
        photo_jacobian_problem_exterior::populate(J, i*2, idx_img, obj_pts, img, cam);
        
        if (prob->no_of_cam_var == 8)
        {
            photo_jacobian_problem_camera_distort::populate(J, i*2, idx_cam, obj_pts, img, cam, prob->img_pts[i]);
        }
        else
        {
            photo_jacobian_problem_camera::populate(J, i*2, idx_cam, obj_pts, img, cam);
        }
        
    }
    
    //print(J);
	return J;
}


/*
 *****************************************************
 * Multi-ray triangulation
 *****************************************************
 */

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
        int img_id = prob->img_pts[i].img_id;
        img img = prob->imgs[img_id];
        
        int cam_id = img.cam_id;
        camera cam = prob->cams[cam_id];  

        int obj_pts_id = prob->img_pts[i].obj_pts_id;
        object_pt obj_pts = prob->obj_pts[obj_pts_id];
        
        point3d pt0;
        pt0.x = x(0);
        pt0.y = x(1);
        pt0.z = x(2);
        
        img_pt pti0;
        backproject(pt0, img, cam, pti0);
        
        r(i*2)   = pti0.x - prob->img_pts[i].x;
        r(i*2+1) = pti0.y - prob->img_pts[i].y;
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
        int img_id = prob->img_pts[i].img_id;
        img img = prob->imgs[img_id];
        
        int cam_id = img.cam_id;
        camera cam = prob->cams[cam_id];        
        
        photo_jacobian_problem_triang::populate(J, i*2, 0, x, img,  cam);
    }
    
    //print(J);
	return J;
}

/*
 *****************************************************
 * Backprojection
 *****************************************************
 */


int backproject(const object_pt &obj_pt, const img &img, const camera &cam, img_pt &pti)
{
    point3d pt3d;
    pt3d.x = obj_pt.x;
    pt3d.y = obj_pt.y;
    pt3d.z = obj_pt.z;    
    return backproject(pt3d, img, cam, pti);
}

int backproject(const point3d &pt, const img &img, const camera &cam, img_pt &pti)
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
    pti.img_id = img.id;
    
    if (cam.cam_type == CAM_TYPE_DISTORTED)
    {   
        double x_hat = pti.x - cam.cx;
        double y_hat = pti.y - cam.cy;
        double r = sqrt(pow(x_hat, 2) + pow(y_hat, 2));
        pti.x -= cam.k1 * pow(r, 3) + cam.k2 * pow(r, 5) + cam.k3 * pow(r, 7) + cam.p1*(pow(r,2) + 2*pow(x_hat,2)) + 2*cam.p2*x_hat*y_hat;
        pti.y -= cam.k1 * pow(r, 3) + cam.k2 * pow(r, 5) + cam.k3 * pow(r, 7) + 2*cam.p1*x_hat*y_hat + cam.p2*(pow(r,2) + 2*pow(y_hat,2)) ;
    }
    
    return 0;
}

int init_problem(problem &prob)
{
    prob.n_unknown_cams = 0;
    for (size_t i = 0; i < prob.cams.size(); i++)
    {
        if (prob.cams[i].type == var_type::UNKNOWN) prob.n_unknown_cams++;
    }   
    
    prob.n_unknown_obj_pts = 0;
    for (size_t i = 0; i < prob.obj_pts.size(); i++)
    {
        if (prob.obj_pts[i].type == var_type::UNKNOWN) prob.n_unknown_obj_pts++;
    }   

    prob.n_unknown_img_pts = 0;
    for (size_t i = 0; i < prob.img_pts.size(); i++)
    {
        if (prob.img_pts[i].type == var_type::UNKNOWN) prob.n_unknown_img_pts++;
    }   
    
    prob.n_unknown_imgs = 0;
    for (size_t i = 0; i < prob.imgs.size(); i++)
    {
        if (prob.imgs[i].type == var_type::UNKNOWN) prob.n_unknown_imgs++;
    } 
    
    prob.n_imgs = prob.imgs.size();
    prob.no_of_img_var = 6;
    prob.idx_imgs = 0;

    prob.n_cams = prob.cams.size();    
    
    // TODO: make it flexible
    if (prob.n_cams > 0)
    {
        if (prob.cams[0].cam_type == CAM_TYPE_DISTORTED)
        {
            prob.no_of_cam_var = 8;
        }
        else
        {
            prob.no_of_cam_var = 3;
        }
    }
    
    prob.idx_cams = prob.idx_imgs +  prob.n_imgs*prob.no_of_img_var;

    prob.sum_unknowns = prob.no_of_img_var*prob.n_imgs + prob.n_cams*prob.no_of_cam_var;
    prob.sum_obs = prob.img_pts.size()*2;
            
    return 0;
};


