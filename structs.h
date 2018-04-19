#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>

#define CAM_TYPE_DISTORTED 2

typedef enum 
{
    KNOWN = 1,
    UNKNOWN = 2
} var_type;

struct point2d
{
    double x;
    double y;
};

struct point3d
{
    double x;
    double y;
    double z;
};

struct object_pt
{
    int id;
    int tid; 
    size_t xid;
    
    double x;
    double y;
    double z;
    var_type type;
};

struct orient3d
{
    double omega;
    double phi;
    double kappa;
};

struct camera
{
    int id; 
    int tid; 
    size_t xid;
    
    int cam_type;
    
    double f;
    double cx;
    double cy; 

    // C++11
    double dc = 0;
    double k1 = 0;
    double k2 = 0;
    double k3 = 0; 
    double p1 = 0;
    double p2 = 0;

    var_type type;    
   
};

struct img
{
    int id; 
    int tid;  
    size_t xid;
    
    double x;
    double y;
    double z;    
    
    double omega;
    double phi;
    double kappa;

    camera* cam_ptr;
    
    var_type type;
};

struct img_pt
{
    int id;
    int tid; 
    
    double x;
    double y;
        
    img* img_ptr;
    object_pt* obj_pts_ptr;
    
    var_type type;
    
};

struct problem
{
    std::vector<img_pt> img_pts;
    std::vector<img> imgs;
    std::vector<object_pt> obj_pts;
    std::vector<camera> cams;
    
    size_t n_unknown_cams;
    size_t n_unknown_obj_pts; 
    size_t n_unknown_img_pts;
    size_t n_unknown_imgs;
    
    size_t sum_unknowns;
    size_t sum_obs;

    size_t start_idx_imgs;
    size_t end_idx_imgs;

    size_t start_idx_cams;
    size_t end_idx_cams;

    size_t start_idx_obj_pts;
    size_t end_idx_obj_pts;
};

#endif
