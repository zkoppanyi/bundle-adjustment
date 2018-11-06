#ifndef CORE_H
#define CORE_H

#include <Eigen/Dense>
#include <vector>
#include "structs.h"
#include "optim.h"

#define MAX_MATRIX_SIZE (100000000/sizeof(double)) // 0.1 GB

int init_problem(problem &prob);

int backproject(object_pt &obj_pt, img &img, const camera &cam, img_pt &pti);
int backproject(const point3d &pt, img &img, const camera &cam, img_pt &pti);
int backproject(object_pt &obj_pt, img &img, img_pt &pti);
int backproject(const point3d &pt, img &img, img_pt &pti);

struct problem_result
{
    optimizer_result optimizer_result;
    double residual;
};
int multiray_triangulate(const problem &prob, point3d &sol, problem_result &result);

struct stochastic_params
{
    int status;
    char status_msg[250];
    
    double sigma_0;
    Eigen::MatrixXd Mxx;
    Eigen::MatrixXd Mll;
};

void calc_stochastic(const optimizer_result &result, stochastic_params &params);

int bundle_adjustment(problem &prob, problem_result &result);
Eigen::VectorXd bundle_adjustment_fn(Eigen::VectorXd x, void* params);
int bundle_adjustment_jacobian(Eigen::VectorXd x, Eigen::SparseMatrix<double>& J, void* params);

#endif