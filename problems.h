#ifndef PROBLEMS_H
#define PROBLEMS_H

#include "matlab.h"
#include "core.h"

template<bool _cX, bool _cY, bool _cZ, bool _cX0, bool _cY0, bool _cZ0, bool _comega, bool _cphi, bool _ckappa, bool _cf, bool _ccx, bool _ccy, bool _ck1, bool _ck2, bool _ck3, bool _cp1, bool _cp2>
struct photo_jacobian_problem 
{   
    static const bool cX = _cX;
    static const bool cY = _cY;
    static const bool cZ = _cZ;
    
    static const bool cX0 = _cX0;
    static const bool cY0 = _cY0;
    static const bool cZ0 = _cZ0;
    static const bool comega = _comega;
    static const bool cphi = _cphi;
    static const bool ckappa = _ckappa;
    
    static const bool cf = _cf;
    static const bool ccx = _ccx;
    static const bool ccy = _ccy;

    static const bool ck1 = _ck1;
    static const bool ck2 = _ck2;
    static const bool ck3 = _ck3;
    static const bool cp1 = _cp1;
    static const bool cp2 = _cp2;
    
    
    // Populater with Eigen::VectorXd
    void static populate(Eigen::MatrixXd &A, size_t i, size_t k, const Eigen::VectorXd &pt0, const img& img0, const camera& cam0)
    {
        if ((cf == true) ||  (ccx == true) || (ccy == true))
        {
            // you cannot call this function on distorted camera model, because it requires to know the image point;
            // use the other populate method
            ASSERT(cam0.cam_type != CAM_TYPE_DISTORTED); 
        } 
        
        populate(A, i, k, pt0(0), pt0(1), pt0(2), img0.x, img0.y, img0.z, img0.omega, img0.phi, img0.kappa, cam0.f, cam0.cx, cam0.cy);
    }

    // Populater with Eigen::VectorXd and image point
    void static populate(Eigen::MatrixXd &A, size_t i, size_t k, const Eigen::VectorXd &pt0, const img& img0, const camera& cam0, const img_pt &pti0)
    {
        populate(A, i, k, pt0(0), pt0(1), pt0(2), img0.x, img0.y, img0.z, img0.omega, img0.phi, img0.kappa, cam0.f, cam0.cx, cam0.cy, pti0.x, pti0.y, cam0.k1, cam0.k2, cam0.k3, cam0.p1, cam0.p2);
    }

    // Populater with point3d
    void static populate(Eigen::MatrixXd &A, size_t i, size_t k, const point3d &pt0, const img &img0, const camera &cam0)
    {
        if ((cf == true) ||  (ccx == true) || (ccy == true))
        {
            // you cannot call this function on distorted camera model, because it requires to know the image point;
            // use the other populate method
            ASSERT(cam0.cam_type != CAM_TYPE_DISTORTED); 
        }
        
        populate(A, i, k, pt0.x, pt0.y, pt0.z, img0.x, img0.y, img0.z, img0.omega, img0.phi, img0.kappa, cam0.f, cam0.cx, cam0.cy);
    }
    
    // Populater with point3d and image point
    void static populate(Eigen::MatrixXd &A, size_t i, size_t k, const point3d &pt0, const img &img0, const camera &cam0, const img_pt &pti0)
    {
        if ((cf == true) ||  (ccx == true) || (ccy == true))
        {
            // you cannot call this function on distorted camera model, because it requires to know the image point;
            // use the other populate method
            ASSERT(cam0.cam_type != CAM_TYPE_DISTORTED); 
        }
                
        populate(A, i, k, pt0.x, pt0.y, pt0.z, img0.x, img0.y, img0.z, img0.omega, img0.phi, img0.kappa, cam0.f, cam0.cx, cam0.cy, pti0.x, pti0.y, cam0.k1, cam0.k2, cam0.k3, cam0.p1, cam0.p2);
    }
    
    // Populater with object_pt
    void static populate(Eigen::MatrixXd &A, size_t i, size_t k, const object_pt &pt0, const img &img0, const camera &cam0)
    {
        if ((cf == true) ||  (ccx == true) || (ccy == true))
        {
            // you cannot call this function on distorted camera model, because it requires to know the image point;
            // use the other populate method
            ASSERT(cam0.cam_type != CAM_TYPE_DISTORTED); 
        }
        
        populate(A, i, k, pt0.x, pt0.y, pt0.z, img0.x, img0.y, img0.z, img0.omega, img0.phi, img0.kappa, cam0.f, cam0.cx, cam0.cy);
    }

    // Populater with object_pt and image point
    void static populate(Eigen::MatrixXd &A, size_t i, size_t k, const object_pt &pt0, const img &img0, const camera &cam0, const img_pt &pti0)
    {
        populate(A, i, k, pt0.x, pt0.y, pt0.z, img0.x, img0.y, img0.z, img0.omega, img0.phi, img0.kappa, cam0.f, cam0.cx, cam0.cy, pti0.x, pti0.y, cam0.k1, cam0.k2, cam0.k3, cam0.p1, cam0.p2);
    }


    void static populate(Eigen::MatrixXd &A, size_t i, size_t k_in, 
                            double X, double Y, double Z, double X0, double Y0, double Z0, double omega, double phi, double kappa, 
                            double f, double cx, double cy, double x = 0, double y = 0, double k1 = 0, double k2 = 0, double k3 = 0, double p1 = 0, double p2 = 0)
    {

        // X component
        int k = (int)k_in-1;
        
        //printf("%.3f, %.3f, %.3f, %.3f, %.3f\n", x, y, k1, k2, k3);
        bool is_dist_params = false;
        double r2 = pow(cx - x, 2) + pow(cy - y,2);
        double r = sqrt(r2);
        if ( (ck1 == true) || (ck2 == true) || (ck3 == true) || (cp1 == true) || (cp2 == true) )
        {           
            is_dist_params = true;
        }
        
        // dX
        if(cX == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*(sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*cos(kappa)*cos(phi))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // dY
        if(cY == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*cos(phi)*sin(kappa))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) + (f*(cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // dZ
        if(cZ == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*cos(omega)*cos(phi)*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*sin(phi))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // dX0
        if(cX0 == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*cos(kappa)*cos(phi))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*(sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // dY0
        if(cY0 == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = - (f*cos(phi)*sin(kappa))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*(cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // dZ0
        if(cZ0 == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*sin(phi))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*cos(omega)*cos(phi)*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // omega
        if(comega == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0))*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // phi
        if(cphi== true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = - (f*(cos(phi)*(Z - Z0) - cos(kappa)*sin(phi)*(X - X0) + sin(kappa)*sin(phi)*(Y - Y0)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0))*(cos(omega)*sin(phi)*(Z - Z0) + cos(kappa)*cos(omega)*cos(phi)*(X - X0) - cos(omega)*cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // kappa
        if(ckappa== true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*(cos(kappa)*cos(phi)*(Y - Y0) + cos(phi)*sin(kappa)*(X - X0)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) + (f*((cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(X - X0) - (sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(Y - Y0))*(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // f
        if(cf == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = -(sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*(X - X0) - cos(phi)*sin(kappa)*(Y - Y0))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // cx
        if(ccx== true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            
            if (is_dist_params)
            {
                A(i,k) = k1*pow(r,2) - 2*p2*(cy - y) + k2*pow(r,4) + k3*pow(r,6) + (cx - x)*(k1*(2*cx - 2*x) + 2*k2*pow(r,2)*(2*cx - 2*x) + 3*k3*pow(r,4)*(2*cx - 2*x)) - p1*(6*cx - 6*x) + 1;
           }
            else
            {
                A(i,k) = 1;
            }
        }

        // cy
        if(ccy == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            
            if (is_dist_params)
            {
                A(i,k) = -2*p2*(cx - x) + (cx - x)*(k1*(2*cy - 2*y) + 2*k2*pow(r, 2)*(2*cy - 2*y) + 3*k3*pow(r, 4)*(2*cy - 2*y)) - p1*(2*cy - 2*y);
            }
            else
            {
                A(i,k) = 0;
            }
        }
        
        // ck1
        if(ck1 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = pow(r,2) * (cx - x);
        }

        // ck2
        if(ck2 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) =  pow(r,4) * (cx - x);
        }
        
        // ck3
        if(ck3 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) =  pow(r,6) * (cx - x);
        }

        // cp1
        if(cp1 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = - 3*pow(cx - x, 2) - pow(cy - y, 2);
        }
        
        // cp2
        if(cp2 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = -2*(cx - x)*(cy - y);
        }

        // Y component
        k = (int)k_in-1;
        i++;

        // X
        if(cX== true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*(sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*(cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // Y
        if(cY == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*(cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*(cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // Z
        if(cZ == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*cos(phi)*sin(omega))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) + (f*cos(omega)*cos(phi)*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // X0
        if(cX0 == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*(cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*(sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // Y0
        if(cY0 == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*(cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*(cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // Z0
        if(cZ0== true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = - (f*cos(phi)*sin(omega))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*cos(omega)*cos(phi)*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // omega
        if(comega == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = f + (f*POW2((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // phi
        if(cphi == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = - (f*(sin(omega)*sin(phi)*(Z - Z0) + cos(kappa)*cos(phi)*sin(omega)*(X - X0) - cos(phi)*sin(kappa)*sin(omega)*(Y - Y0)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0))*(cos(omega)*sin(phi)*(Z - Z0) + cos(kappa)*cos(omega)*cos(phi)*(X - X0) - cos(omega)*cos(phi)*sin(kappa)*(Y - Y0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // kappa
        if(ckappa == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = (f*((cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(X - X0) - (sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(Y - Y0))*((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0)))/POW2((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0)) - (f*((cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(X - X0) - (cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(Y - Y0)))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }
        
        // f
        if(cf == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = -((cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi))*(X - X0) + (cos(kappa)*cos(omega) - sin(kappa)*sin(omega)*sin(phi))*(Y - Y0) - cos(phi)*sin(omega)*(Z - Z0))/((sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi))*(X - X0) + (cos(kappa)*sin(omega) + cos(omega)*sin(kappa)*sin(phi))*(Y - Y0) + cos(omega)*cos(phi)*(Z - Z0));
        }

        // cx
        if(ccx == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );

            if (is_dist_params)
            {
                A(i,k) = (cy - y)*(k1*(2*cx - 2*x) + 2*k2*pow(r,2)*(2*cx - 2*x) + 3*k3*pow(r,4)*(2*cx - 2*x)) - 2*p1*(cy - y) - p2*(2*cx - 2*x);
            }
            else
            {
                A(i,k) = 0;
            }
            
        }

        // cy
        if(ccy == true)
        {
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            
            if (is_dist_params)
            {
                A(i,k) = k1*pow(r,2) - 2*p1*(cx - x) + k2*pow(r,4) + k3*pow(r,6) + (cy - y)*(k1*(2*cy - 2*y) + 2*k2*pow(r,2)*(2*cy - 2*y) + 3*k3*pow(r,4)*(2*cy - 2*y)) - p2*(6*cy - 6*y) + 1; 
            }
            else
            {
                A(i,k) = 1;
            }
        }
        
        // ck1
        if(ck1 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = pow(r,2)*(cy - y);
        }

        // ck2
        if(ck2 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = pow(r,4)*(cy - y);
        }
        
        // ck3
        if(ck3 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = pow(r,6)*(cy - y);
        }

        // cp1
        if(cp1 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = -2*(cx - x)*(cy - y);
        }
        
        // cp2
        if(cp2 == true)
        {        
            k++;
            ASSERT( (0 <= i) && (i < A.rows()) && (0 <= k) && (k < A.cols()) );
            A(i,k) = -1 * pow((cx - x), 2) - 3*pow((cy - y), 2);
        }
        
    }
};

typedef photo_jacobian_problem<true, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false> photo_jacobian_problem_triang;
typedef photo_jacobian_problem<false, false, false, false, false, false, false, false, false, true, true, true, false, false, false, false, false> photo_jacobian_problem_camera;
typedef photo_jacobian_problem<false, false, false, false, false, false, false, false, false, true, true, true, true, true, true, true, true> photo_jacobian_problem_camera_distort;
typedef photo_jacobian_problem<false, false, false, true, true, true, true, true, true, false, false, false, false, false, false, false, false> photo_jacobian_problem_exterior;

        
#endif