#ifndef OPTIM_H
#define OPTIM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

enum stopping_criteria
{
    THRESHOLD_REACHED,
    MAXIMUM_ITERATION_REACHED,
    UNBALANCED_PROBLEM,
    ERROR
};

struct optimizer_result
{
    stopping_criteria stopping_criteria;
    size_t no_of_iterations;
    Eigen::SparseMatrix<double> J;
    Eigen::VectorXd r;
};

Eigen::VectorXd levenberg_marquardt(Eigen::VectorXd (*fn)(Eigen::VectorXd, void*), int(*jacobian)(Eigen::VectorXd,  Eigen::SparseMatrix<double>&, void*), void* params, Eigen::VectorXd x, double TolX, double TolY, optimizer_result &result);
//Eigen::VectorXd simulated_annealing(Eigen::VectorXd (*fn)(Eigen::VectorXd, void*), Eigen::MatrixXd (*jacobian)(Eigen::VectorXd, void*), void* params, Eigen::VectorXd x, const double TolX, const double TolY, optimizer_result &result);
#endif