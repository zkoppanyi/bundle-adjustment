#include "matlab.h"

void print(Eigen::MatrixXd A)
{
    mexPrintf("\nMatrix (%i, %i)\n", A.rows(), A.cols());
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            mexPrintf(" %.5f ", A(i, j));
        }
        mexPrintf("\n");
    }
}

void print(Eigen::VectorXd v)
{
    mexPrintf("\nVector (%i)\n", v.rows());
    for (int i = 0; i < v.rows(); i++)
    {
       mexPrintf(" %.5f \n", v(i));
    }
}

void eigen2mat(const Eigen::MatrixXd &A, mxArray* &ret)
{
    size_t n_rows = A.rows();
    size_t n_cols = A.cols();

    ret = mxCreateDoubleMatrix(n_rows, n_cols, mxREAL);
    double* ret_ptr = mxGetPr(ret);
     
    for(size_t i = 0; i < n_rows; i++)
    {
        for(size_t j = 0; j < n_cols; j++)
        {
            GET(ret_ptr, i, j, n_rows) = A(i, j);
        }
    }   
}

void eigen2mat(double d, mxArray* &ret)
{
    ret = mxCreateDoubleScalar(d);
}

//
// MODE = 1 : read all;
// MODE = 2 : read for backrpojection;
void extract_problem_from_arguments(int nrhs, const mxArray *prhs[], problem &prob, int mode)
{
    int pts_arr_idx = 0;
    int img_arr_idx = 1;
    int obj_pts_arr_idx = 2;
    int camera_arr_idx = 3;
    int param_req = 4;
    
    if (mode == 1)
    {
        // general case
    }
    else if (mode == 2)
    {        
        pts_arr_idx = -1; // ignore
        img_arr_idx = 0;
        obj_pts_arr_idx = 1;
        camera_arr_idx = 2;
        param_req = 3;
    }
    
    if (nrhs < param_req)
    {
        LOG2I("Number of input arguments: ", nrhs);
        LOG2I("Number of required arguments: ", param_req);
        LOG_ERROR("Not enough input parameters!");
        //LOG_ERROR("Not enough input arguments! At least 4 required!\n Parameters: img_pts, imgs, point3ds, cams\n");
        return;
    }
    
    // Load 2d image points
    double* pts_arr = NULL;
    size_t pts_m = 0;
    size_t pts_n = 0;
    if (pts_arr_idx >= 0)
    {
        pts_arr = mxGetPr(prhs[pts_arr_idx]);
        pts_m = mxGetN(prhs[pts_arr_idx]);
        pts_n = mxGetM(prhs[pts_arr_idx]);

        if (pts_m < 6)
        {
            LOG_ERROR("First argument has to be a matrix of 2D image points: [nx6] (one row: ID, x, y, cam_id, obj_pts_id, type) \n");
            return;
        }
    }
   
    //LOG2I("No. of points: ", pts_n);  
    
    // Load images    
    double* img_arr = NULL;
    size_t img_m = 0;
    size_t img_n = 0;
    if (img_arr_idx >= 0)
    {
        img_arr = mxGetPr(prhs[img_arr_idx]);
        img_m = mxGetN(prhs[img_arr_idx]);
        img_n = mxGetM(prhs[img_arr_idx]);

        if (img_m < 9)
        {
            LOG_ERROR("Second argument has to be the images with at least 1+6+2 columns (one row: ID, X0, Y0, Z0, Omega, Phi, Kappa, cam_id, type)\n");
            return;
        }
    }

    // Load 3d points
    double* obj_pts_arr = NULL;
    size_t obj_pts_m = 0;
    size_t obj_pts_n = 0;
    if (obj_pts_arr_idx >= 0)
    {            
        obj_pts_arr = mxGetPr(prhs[obj_pts_arr_idx]);
        obj_pts_m = mxGetN(prhs[obj_pts_arr_idx]);
        obj_pts_n = mxGetM(prhs[obj_pts_arr_idx]);

        if (obj_pts_m < 4)
        {
            LOG_ERROR("Third argument has to be the 3d points with at least 5 columns (one row: ID, X, Y, Z, type)\n");
            return;
        }
    }
    
    // Load cameras    
    double* camera_arr = NULL;    
    size_t camera_m = 0;
    size_t camera_n = 0;
    if (camera_arr_idx >= 0)
    {  
        camera_arr = mxGetPr(prhs[camera_arr_idx]);
        camera_m = mxGetN(prhs[camera_arr_idx]);
        camera_n = mxGetM(prhs[camera_arr_idx]);

        if (camera_m < 6)
        {
            LOG_ERROR("Fourth argument has to be the camera's matrix with at least 5 columns (one row: ID, cam_type, f, cx, cy, [k1, k2, k3, p1, p2] type)\n");
            return;
        }
    }
    
    // create problem   
    
    // read up cameras
    for (size_t i = 0; i < camera_n; i++)
    {
        camera cam;
        cam.id        = (int)GET(camera_arr, i, 0, camera_n) - 1;
        cam.cam_type  = (int)GET(camera_arr, i, 1, camera_n);
        cam.f         =      GET(camera_arr, i, 2, camera_n);
        cam.cx        =      GET(camera_arr, i, 3, camera_n);
        cam.cy        =      GET(camera_arr, i, 4, camera_n);
        
        size_t k = 5;
        if (cam.cam_type == CAM_TYPE_DISTORTED)
        {
            if (camera_m < 11)
            {
                LOG_ERROR("cam_type indicates using distortion models but they are not defined in cams!"); 
                return;
            }
        }
        
        if (camera_m >= 10)
        {            
            cam.k1  =      GET(camera_arr, i, 5, camera_n);
            cam.k2  =      GET(camera_arr, i, 6, camera_n);
            cam.k3  =      GET(camera_arr, i, 7, camera_n);
            cam.p1  =      GET(camera_arr, i, 8, camera_n);
            cam.p2  =      GET(camera_arr, i, 9, camera_n);                
            k = 10;            
        }

        
        int type_id   = (int)GET(camera_arr, i, k, camera_n);

        switch(type_id)
        {
            case static_cast<int>(var_type::KNOWN)   : cam.type = var_type::KNOWN; break;
            case static_cast<int>(var_type::UNKNOWN) : cam.type = var_type::UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in cams!"); 
                return;
                break;
        }
        
        if ((cam.id < 0) || cam.id >= camera_n)
        {
            mexPrintf("Invalid ID: %i at (%i, %i) in cams\n", cam.id, (int)i+1, 5); 
            LOG_ERROR("Invalid camera ID in cams!"); 
            return;
        }
        
        prob.cams.push_back(cam);
    }
    
    // read up object space points    
    for (size_t i = 0; i < obj_pts_n; i++)
    {       
        object_pt pt3d;
        pt3d.id      = (int)GET(obj_pts_arr, i, 0, obj_pts_n) - 1;
        pt3d.x       = GET(obj_pts_arr, i, 1, obj_pts_n);
        pt3d.y       = GET(obj_pts_arr, i, 2, obj_pts_n);
        pt3d.z       = GET(obj_pts_arr, i, 3, obj_pts_n);
        int type_id  = (int)GET(obj_pts_arr, i, 4, obj_pts_n);    
        
         switch(type_id)
        {
            case static_cast<int>(var_type::KNOWN)   : pt3d.type = var_type::KNOWN; break;
            case static_cast<int>(var_type::UNKNOWN) : pt3d.type = var_type::UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in obj_pts!"); 
                return;
                break; 
        }
         
        if ((pt3d.id < 0) || pt3d.id >= obj_pts_n)
        {
            mexPrintf("Invalid id: %i at (%i, %i) in pt3ds\n", pt3d.id, i+1, 4);
            LOG_ERROR("Invalid object point ID in pt3ds!"); 
            return;
        }
        
        prob.obj_pts.push_back(pt3d);
    }
    
    // read up images    
    for (size_t i = 0; i < img_n; i++)
    {       
        img img;
        img.id         = (int)GET(img_arr, i, 0, img_n) - 1;
        img.x          =      GET(img_arr, i, 1, img_n);
        img.y          =      GET(img_arr, i, 2, img_n);
        img.z          =      GET(img_arr, i, 3, img_n);
        img.omega      =      GET(img_arr, i, 4, img_n);
        img.phi        =      GET(img_arr, i, 5, img_n);
        img.kappa      =      GET(img_arr, i, 6, img_n);
        img.cam_id     = (int)GET(img_arr, i, 7, img_n) - 1;
        int type_id    = (int)GET(img_arr, i, 8, img_n);
        
        switch(type_id)
        {
            case static_cast<int>(var_type::KNOWN)   : img.type = var_type::KNOWN; break;
            case static_cast<int>(var_type::UNKNOWN) : img.type = var_type::UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in imgs!"); 
                return;
                break;        
        }
        
        if ((img.id < 0) || img.id >= img_n)
        {
            mexPrintf("Invalid id: %i at (%i, %i) in imgs\n", img.id, i+1, 4);
            LOG_ERROR("Invalid object point ID in imgs!"); 
            return;
        }
        
        if ((img.cam_id < 0) || img.cam_id >= camera_n)
        {
            mexPrintf("Invalid cam_id: %i at (%i, %i) in imgs\n", img.id, i+1, 4);
            LOG_ERROR("Invalid camera ID in imgs!"); 
            return;
        }
        
        prob.imgs.push_back(img);
    }
    
    // read up image points    
    for (size_t i = 0; i < pts_n; i++)
    {       
        img_pt pti;
        pti.id         = (int)GET(pts_arr, i, 0, pts_n) - 1;
        pti.x          =      GET(pts_arr, i, 1, pts_n);
        pti.y          =      GET(pts_arr, i, 2, pts_n);
        pti.img_id     = (int)GET(pts_arr, i, 3, pts_n) - 1;
        pti.obj_pts_id = (int)GET(pts_arr, i, 4, pts_n) - 1;
        int type_id    = (int)GET(pts_arr, i, 5, pts_n);
        
        switch(type_id)
        {
            case static_cast<int>(var_type::KNOWN)   : pti.type = var_type::KNOWN; break;
            case static_cast<int>(var_type::UNKNOWN) : pti.type = var_type::UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in img_pts!"); 
                return;
                break;        
        }
        
        if ((pti.img_id < 0) || pti.img_id >= img_n)
        {
            mexPrintf("Invalid img_id: %i at (%i, %i) in img_id\n", pti.id, i+1, 4);
            LOG_ERROR("Invalid camera ID in img_pts!"); 
            return;
        }

        if ((pti.obj_pts_id < 0) || pti.obj_pts_id >= obj_pts_n)
        {
            mexPrintf("Invalid obj_pts_id: %i at (%i, %i) in img_pts\n", pti.obj_pts_id, i+1, 5);
            LOG_ERROR("Invalid object point ID in img_pts!"); 
            return;
        }

        prob.img_pts.push_back(pti);
    }
}

void create_stoch_struct(const optimizer_result &optimizer_result, const stochastic_params &stoch, mxArray* &ret)
{
        mxArray* Mxx;        
        eigen2mat(stoch.Mxx, Mxx); 
        
        mxArray* Mll;        
        eigen2mat(stoch.Mll, Mll); 

        mxArray* J;        
        eigen2mat(optimizer_result.J, J); 

        mxArray* sigma_0;        
        eigen2mat(stoch.sigma_0, sigma_0);

        const char *fieldnames[] = {"Mxx", "Mll", "J", "sigma0"};        
        ret = mxCreateStructMatrix(1,1,4,fieldnames);
        
        mxSetFieldByNumber(ret,0,0, Mxx);
        mxSetFieldByNumber(ret,0,1, Mll);
        mxSetFieldByNumber(ret,0,2, J);
        mxSetFieldByNumber(ret,0,3, sigma_0);   
}

void create_problem_struct(const problem &prob, mxArray* &ret)
{
        // imgs
        size_t n_ret = prob.imgs.size();
        mxArray* img_arr = mxCreateDoubleMatrix(n_ret, 9, mxREAL);
        double* img_ptr = mxGetPr(img_arr);
        for (size_t i = 0; i < n_ret ; i++)
        {
            img img = prob.imgs[i];

            GET(img_ptr, i, 0, n_ret) = img.id + 1;
            GET(img_ptr, i, 1, n_ret) = img.x;
            GET(img_ptr, i, 2, n_ret) = img.y;
            GET(img_ptr, i, 3, n_ret) = img.z;
            GET(img_ptr, i, 4, n_ret) = img.omega;
            GET(img_ptr, i, 5, n_ret) = img.phi;
            GET(img_ptr, i, 6, n_ret) = img.kappa;
            GET(img_ptr, i, 7, n_ret) = img.cam_id + 1;
            GET(img_ptr, i, 8, n_ret) = static_cast<int>(img.type);
                    }   
      
        n_ret = prob.cams.size();
        
        size_t m_cam = 6;        
        bool is_distort = false;
        
        // discover whether we have distorted camera, and return based on that
        for (size_t i = 0; i < n_ret ; i++)
        {
            if (prob.cams[0].cam_type == CAM_TYPE_DISTORTED) 
            {
                is_distort = true;
                m_cam = 11;
                break;
            }
        }
                
        mxArray* cams_arr = mxCreateDoubleMatrix(n_ret, m_cam, mxREAL);
        double* cams_ptr = mxGetPr(cams_arr);
        for (size_t i = 0; i < n_ret ; i++)
        {
            camera cam = prob.cams[i];
            GET(cams_ptr, i, 0, n_ret) = cam.id + 1;
            GET(cams_ptr, i, 1, n_ret) = cam.type;
            GET(cams_ptr, i, 2, n_ret) = cam.f;
            GET(cams_ptr, i, 3, n_ret) = cam.cx;
            GET(cams_ptr, i, 4, n_ret) = cam.cy;
            
            if (is_distort)
            {
                GET(cams_ptr, i, 5, n_ret) = cam.k1;
                GET(cams_ptr, i, 6, n_ret) = cam.k2;
                GET(cams_ptr, i, 7, n_ret) = cam.k3;
                GET(cams_ptr, i, 8, n_ret) = cam.p1;
                GET(cams_ptr, i, 9, n_ret) = cam.p2;
                GET(cams_ptr, i, 10, n_ret) = static_cast<int>(cam.type);
            }
            else
            {
                GET(cams_ptr, i, 5, n_ret) = static_cast<int>(cam.type);
            }
        }   
                
        n_ret = prob.img_pts.size();
        mxArray* img_pts_arr = mxCreateDoubleMatrix(n_ret, 6, mxREAL);
        double* img_pts_ptr = mxGetPr(img_pts_arr);
        for (size_t i = 0; i < n_ret ; i++)
        {
            img_pt pt = prob.img_pts[i];
            GET(img_pts_ptr, i, 0, n_ret) = pt.id + 1;
            GET(img_pts_ptr, i, 1, n_ret) = pt.x;
            GET(img_pts_ptr, i, 2, n_ret) = pt.y;
            GET(img_pts_ptr, i, 3, n_ret) = pt.img_id + 1;
            GET(img_pts_ptr, i, 4, n_ret) = pt.obj_pts_id + 1;
            GET(img_pts_ptr, i, 5, n_ret) = static_cast<int>(pt.type);
        }   

        n_ret = prob.obj_pts.size();
        mxArray* obj_pts_arr = mxCreateDoubleMatrix(n_ret, 5, mxREAL);
        double* obj_pts_ptr = mxGetPr(obj_pts_arr);
        for (size_t i = 0; i < n_ret ; i++)
        {
            object_pt pt = prob.obj_pts[i];

            GET(obj_pts_ptr, i, 0, n_ret) = pt.id + 1;
            GET(obj_pts_ptr, i, 1, n_ret) = pt.x;
            GET(obj_pts_ptr, i, 2, n_ret) = pt.y;
            GET(obj_pts_ptr, i, 3, n_ret) = pt.z;
            GET(obj_pts_ptr, i, 4, n_ret) = static_cast<int>(pt.type);
        }   

        const char *fieldnames[] = {"img_pts", "imgs", "obj_pts", "cams"};        
        ret = mxCreateStructMatrix(1,1,4,fieldnames);
        
        mxSetFieldByNumber(ret,0,0, img_pts_arr);
        mxSetFieldByNumber(ret,0,1, img_arr);
        mxSetFieldByNumber(ret,0,2, obj_pts_arr);
        mxSetFieldByNumber(ret,0,3, cams_arr);
}
