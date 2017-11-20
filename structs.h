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
    
    double x;
    double y;
    double z;    
    
    double omega;
    double phi;
    double kappa;

    int cam_id;
    var_type type;
};

struct img_pt
{
    int id;
    
    double x;
    double y;
    
    int img_id;
    int obj_pts_id;
    
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

    size_t n_imgs;
    size_t no_of_img_var;
    size_t idx_imgs;

    size_t n_cams;    
    size_t no_of_cam_var;
    size_t idx_cams;
};

#endif
