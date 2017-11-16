#include "matlab.h"

#include <vector>

#include "utils.h"
#include "structs.h"
#include "core.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Load problem
    problem prob;    
    extract_problem_from_arguments(nrhs, prhs, prob);
    
    if (nlhs < 1)
    {
        LOG_ERROR("Need a return value!\n");
        return;
    }
    
    /*init_problem(prob);    
    LOG2I(" Unknown image points: ", prob.n_unknown_img_pts);
    LOG2I("       Unknown images: ", prob.n_unknown_imgs);
    LOG2I("Unknown object points: ", prob.n_unknown_obj_pts);
    LOG2I("      Unknown cameras: ", prob.n_unknown_cams);*/
    
    size_t n_obj_pts = prob.obj_pts.size();
    size_t n_imgs = prob.imgs.size();

    size_t n_ret = n_obj_pts*n_imgs;
    plhs[0] = mxCreateDoubleMatrix(n_ret, 6, mxREAL);
    double* ret = mxGetPr(plhs[0]);
    int k = 0;
    for (size_t i = 0; i < n_obj_pts; i++)
    {
        object_pt obj_pt = prob.obj_pts[i];

        for (size_t j = 0; j < n_imgs; j++)
        {
            img img = prob.imgs[j];

            int cam_id = img.cam_id;
            camera cam = prob.cams[cam_id];
        
            img_pt pti;
            backproject(obj_pt, img, cam, pti);
            pti.id = k;
            pti.obj_pts_id = i;
            pti.type = var_type::KNOWN;

            GET(ret, k, 0, n_ret) = pti.id + 1;
            GET(ret, k, 1, n_ret) = pti.x;
            GET(ret, k, 2, n_ret) = pti.y;
            GET(ret, k, 3, n_ret) = pti.img_id + 1;
            GET(ret, k, 4, n_ret) = pti.obj_pts_id + 1;
            GET(ret, k, 5, n_ret) = static_cast<int>(pti.type);

            /*mexPrintf("x=%.3f y=%.3f z=%.3f\n", pt.x, pt.y, pt.z);
            mexPrintf("x=%.3f, y=%.3f\n", pt2d.x, pt2d.y);*/
            k++;
        }
    }
    
}