#ifndef CORE_H
#define CORE_H

#include <Eigen/Dense>
#include <vector>
#include "structs.h"
#include "optim.h"

int init_problem(problem &prob);

int backproject(const object_pt &obj_pt, const img &img, const camera &cam, img_pt &pti);
int backproject(const point3d &pt, const img &img, const camera &cam, img_pt &pti);

struct problem_result
{
    optimizer_result optimizer_result;
    double residual;
};
int multiray_triangulate(const problem &prob, point3d &sol, problem_result &result);

struct stochastic_params
{
    double sigma_0;
    Eigen::MatrixXd Mxx;
    Eigen::MatrixXd Mll;
};
void calc_stochastic(const optimizer_result &result, stochastic_params &params);

int bundle_adjustment(const problem &prob, problem &sol, problem_result &result);

#endif