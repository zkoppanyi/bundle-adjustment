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
// MODE = 2 : read for backprojection
int extract_problem_from_arguments(int nrhs, const mxArray *prhs[], problem &prob, int mode)
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
        return -1;
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
            return -1;
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
            return -1;
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
            return -1;
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
            return -1;
        }
    }
    
    // create problem   
    
    // read up cameras
    int ctid = 0;
    for (size_t i = 0; i < camera_n; i++)
    {
        camera cam;
        cam.id        = (int)GET(camera_arr, i, 0, camera_n);
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
                return -1;
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
            case static_cast<int>(KNOWN)   : cam.type = KNOWN; break;
            case static_cast<int>(UNKNOWN) : cam.type = UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in cams!"); 
                return -1;
        }
                
        cam.tid = ctid;
        prob.cams.push_back(cam);
        ctid++;
    }
    
    // read up object space points    
    int otid = 0;
    for (size_t i = 0; i < obj_pts_n; i++)
    {       
        object_pt pt3d;
        pt3d.id      = (int)GET(obj_pts_arr, i, 0, obj_pts_n);
        pt3d.x       = GET(obj_pts_arr, i, 1, obj_pts_n);
        pt3d.y       = GET(obj_pts_arr, i, 2, obj_pts_n);
        pt3d.z       = GET(obj_pts_arr, i, 3, obj_pts_n);
        int type_id  = (int)GET(obj_pts_arr, i, 4, obj_pts_n);    
        
        switch(type_id)
        {
            case static_cast<int>(KNOWN)   : pt3d.type = KNOWN; break;
            case static_cast<int>(UNKNOWN) : pt3d.type = UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in obj_pts!"); 
                return -1;
        }
                         
        pt3d.tid = otid;
        prob.obj_pts.push_back(pt3d);
        otid++;
    }
    
    // read up images 
    int itid = 0;
    for (size_t i = 0; i < img_n; i++)
    {       
        img img;
        img.id         = (int)GET(img_arr, i, 0, img_n);
        img.x          =      GET(img_arr, i, 1, img_n);
        img.y          =      GET(img_arr, i, 2, img_n);
        img.z          =      GET(img_arr, i, 3, img_n);
        img.omega      =      GET(img_arr, i, 4, img_n);
        img.phi        =      GET(img_arr, i, 5, img_n);
        img.kappa      =      GET(img_arr, i, 6, img_n);
        
         // finding camera id
        int cam_id     = (int)GET(img_arr, i, 7, img_n);
        
        camera* fcam = NULL;
        for(std::vector<camera>::size_type k = 0; k != prob.cams.size(); k++) 
        {
             if (prob.cams[k].id == cam_id) fcam = &prob.cams[k];
        }        
        img.cam_ptr = fcam;
        
        if (img.cam_ptr == NULL)
        {
            mexPrintf("Invalid cam_id: %i at (%i, %i) in imgs\n", cam_id, i+1, 4);
            LOG_ERROR("Invalid camera ID in imgs!"); 
            return -1;
        }
        
         // finding type
        int type_id    = (int)GET(img_arr, i, 8, img_n);
        
        switch(type_id)
        {
            case static_cast<int>(KNOWN)   : img.type = KNOWN; break;
            case static_cast<int>(UNKNOWN) : img.type = UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in imgs!"); 
                return -1;
        }              
        
        img.tid = itid;
        prob.imgs.push_back(img);
        itid++;
    }

    // read up image points    
    int ptid = 0;
    for (size_t i = 0; i < pts_n; i++)
    {       
        img_pt pti;
        pti.id         = (int)GET(pts_arr, i, 0, pts_n);
        pti.x          =      GET(pts_arr, i, 1, pts_n);
        pti.y          =      GET(pts_arr, i, 2, pts_n);
        
        // finding image id
        int img_id     = (int)GET(pts_arr, i, 3, pts_n);
        
        img* fimg = NULL;
        for(std::vector<img>::size_type k = 0; k != prob.imgs.size(); k++) 
        {
             if (prob.imgs[k].id == img_id) 
             {
                 fimg = &(prob.imgs[k]);
                 break;
             }
        }        
        pti.img_ptr = fimg;
        
        if (pti.img_ptr == NULL)
        {
            mexPrintf("Invalid img_id=%i at (%i, %i) in img_pts!\n", img_id, i+1, 4);
            LOG_ERROR("Invalid image ID in img_pts!"); 
            return -1;
        }
        
        // finding object point id
        int obj_pts_id = (int)GET(pts_arr, i, 4, pts_n);
        
        object_pt* fobj_pt_ptr = NULL;
        for(std::vector<object_pt>::size_type k = 0; k != prob.obj_pts.size(); k++) 
        {
             if (prob.obj_pts[k].id == obj_pts_id) 
             {
                 fobj_pt_ptr = &prob.obj_pts[k];
                 break;
             }
        }        
        pti.obj_pts_ptr = fobj_pt_ptr;
        
        if (pti.obj_pts_ptr == NULL)
        {
            mexPrintf("Invalid obj_pts_id: %i at (%i, %i) in img_pts\n", obj_pts_id, i+1, 5);
            LOG_ERROR("Invalid object point ID in img_pts!"); 
            return -1;
        }
        
        // finding type
        int type_id    = (int)GET(pts_arr, i, 5, pts_n);
        
        switch(type_id)
        {
            case static_cast<int>(KNOWN)   : pti.type = KNOWN; break;
            case static_cast<int>(UNKNOWN) : pti.type = UNKNOWN; break;
            default: 
                LOG2I("Invalid record type at ", (int)i+1); 
                LOG2I("Input record type: ", type_id); 
                LOG_ERROR("Invalid record type in img_pts!"); 
                return -1;
        }
        
        pti.tid = ptid;
        prob.img_pts.push_back(pti);
        ptid++;
    }
    
    return 0;
}

void create_stoch_struct(const optimizer_result &optimizer_result, const stochastic_params &stoch, mxArray* &ret)
{
        if (stoch.status == 0)
        {        
            /*mxArray* Mxx;        
            eigen2mat(stoch.Mxx, Mxx); 

            mxArray* Mll;        
            eigen2mat(stoch.Mll, Mll); 

            mxArray* J;        
            eigen2mat(optimizer_result.J, J); 

            mxArray* sigma_0;        
            eigen2mat(stoch.sigma_0, sigma_0);

            const char *fieldnames[] = {"Mxx", "Mll", "J", "sigma0", "status", "status_msg"};        
            ret = mxCreateStructMatrix(1,1,4,fieldnames);

            mxSetFieldByNumber(ret,0,0, Mxx);
            mxSetFieldByNumber(ret,0,1, Mll);
            mxSetFieldByNumber(ret,0,2, J);
            mxSetFieldByNumber(ret,0,3, sigma_0); 
            mxSetFieldByNumber(ret,0,4, mxCreateDoubleScalar(0)); 
            mxSetFieldByNumber(ret,0,5, mxCreateString("OK")); */

            mxArray* J;        
            eigen2mat(optimizer_result.J, J); 
            const char *fieldnames[] = {"J", "status", "status_msg"};        
            ret = mxCreateStructMatrix(1,1,3,fieldnames);
            mxSetFieldByNumber(ret, 0,0, J);
            mxSetFieldByNumber(ret, 0,1, mxCreateDoubleScalar(-1));
            mxSetFieldByNumber(ret, 0,2, mxCreateString(stoch.status_msg));
        }
        else
        {
            const char *fieldnames[] = {"status", "status_msg"};        
            ret = mxCreateStructMatrix(1,1,2,fieldnames);
            mxSetFieldByNumber(ret, 0,0, mxCreateDoubleScalar(-1));
            mxSetFieldByNumber(ret, 0,1, mxCreateString(stoch.status_msg));            
        }
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

            GET(img_ptr, i, 0, n_ret) = img.id;
            GET(img_ptr, i, 1, n_ret) = img.x;
            GET(img_ptr, i, 2, n_ret) = img.y;
            GET(img_ptr, i, 3, n_ret) = img.z;
            GET(img_ptr, i, 4, n_ret) = img.omega;
            GET(img_ptr, i, 5, n_ret) = img.phi;
            GET(img_ptr, i, 6, n_ret) = img.kappa;
            GET(img_ptr, i, 7, n_ret) = img.cam_ptr->id;
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
                
        // cams
        mxArray* cams_arr = mxCreateDoubleMatrix(n_ret, m_cam, mxREAL);
        double* cams_ptr = mxGetPr(cams_arr);
        for (size_t i = 0; i < n_ret ; i++)
        {
            camera cam = prob.cams[i];
            GET(cams_ptr, i, 0, n_ret) = cam.id;
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
                
        //img_pts
        n_ret = prob.img_pts.size();
        mxArray* img_pts_arr = mxCreateDoubleMatrix(n_ret, 6, mxREAL);
        double* img_pts_ptr = mxGetPr(img_pts_arr);
        for (size_t i = 0; i < n_ret ; i++)
        {
            img_pt pt = prob.img_pts[i];
            GET(img_pts_ptr, i, 0, n_ret) = pt.id;
            GET(img_pts_ptr, i, 1, n_ret) = pt.x;
            GET(img_pts_ptr, i, 2, n_ret) = pt.y;
            GET(img_pts_ptr, i, 3, n_ret) = pt.img_ptr->id;
            GET(img_pts_ptr, i, 4, n_ret) = pt.obj_pts_ptr->id;
            GET(img_pts_ptr, i, 5, n_ret) = static_cast<int>(pt.type);
        }   

        n_ret = prob.obj_pts.size();
        mxArray* obj_pts_arr = mxCreateDoubleMatrix(n_ret, 5, mxREAL);
        double* obj_pts_ptr = mxGetPr(obj_pts_arr);
        for (size_t i = 0; i < n_ret ; i++)
        {
            object_pt pt = prob.obj_pts[i];

            GET(obj_pts_ptr, i, 0, n_ret) = pt.id;
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
