#include <Eigen/IterativeLinearSolvers>
        
#include "optim.h"
#include <cmath>
#include <iostream>
//#include <random>
//#include "matlab.h"

using namespace std;
using namespace Eigen;

bool any_greater_than(VectorXd v1, VectorXd v2)
{
    for (int ik = 0; ik < v1.size(); ik++)
    {
        if (fabs(v1(ik)) >= v2(ik))
        {
            return true;
        }
    }

    return false;
}

/*
 *****************************************************
 * Levenberg-Marquardt optimizer with known Jacobians
 *****************************************************
 */

VectorXd levenberg_marquardt(VectorXd (*fn)(VectorXd, void*), MatrixXd (*jacobian)(VectorXd, void*), void* params, VectorXd x, const double TolX, const double TolY, optimizer_result &result)
{
	const unsigned int MaximumIterationNumber = 500;

	// Residual at starting point
    VectorXd r = fn(x, params);
    double S = r.dot(r);

    size_t lr = r.size();
    size_t lx = x.size();

    double jepsx = TolX;

    MatrixXd J = jacobian(x, params);
    //Debug.WriteLine("Jacobian: " + J);

    int nfJ = 2;
    MatrixXd A = J.transpose() * J;
    VectorXd v = J.transpose() * r;

    // Automatic scaling
    MatrixXd D = A.diagonal().asDiagonal();
    size_t Ddim = min(D.rows(), D.cols());
    for (int i = 0; i < Ddim; i++)
    {
        if (D(i, i) == 0) D(i, i) = 1.0;
    }

    double Rlo = 0.25;
    double Rhi = 0.75;
    double l = 1;
    double lc = 0.75;
    unsigned int cnt = 0;

    VectorXd epsx = VectorXd::Ones(lx) * TolX;
    VectorXd epsy = VectorXd::Ones(lr) * TolY;

    VectorXd d = VectorXd::Ones(lx) * TolX;

    //while (cnt < MaximumIterationNumber)
    while ((cnt < MaximumIterationNumber) && any_greater_than(d, epsx) && (any_greater_than(r, epsy)))
    {
        // negative solution increment
        //d = SolveLinearEquationSystem(A.Add((l.Multiply(D))), v);
        d = (A + l*D).ldlt().solve(v); // use ldlt for stability.
        //d = (A + l*D).llt().solve(v); 
        
        //d = (A + l*D).householderQr().solve(v);
        
        /*LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
        lscg.compute( (A + l*D).sparseView() );
        d = lscg.solve(v);*/
        
        VectorXd xd = x - d;
        VectorXd rd = fn(xd, params);

        nfJ = nfJ + 1;
        double Sd = rd.dot(rd);
        double dS = d.dot(2*v - A*d); // predicted reduction

        double R = (S - Sd) / dS;
        if (R > Rhi)
        {
            l = l / 2; // halve lambda if R too high
            if (l < lc) l = 0;
        }
        else if (R < Rlo)  // find new nu if R too low
        {
            double nu = (Sd - S) / (d.dot(v)) + 2;
            if (nu < 2)
                nu = 2;
            else if (nu > 10)
            {
                nu = 10;
            }

            if (l == 0)
            {
                VectorXd diag = A.inverse().diagonal();
                
                /*double max_pos = diag.Max();
                double max_neg = Math.Abs(diag.Min());*/
                double max_pos = diag(0);
                double max_neg = diag(0);
                for (int k=0; k<diag.size();k++)
                {
                	max_pos = max_pos > diag(k) ? max_pos : diag(k);
                	max_neg = max_neg < diag(k) ? max_neg : diag(k);
                }
				max_neg = abs(max_neg);

                double abs_max = max_pos > max_neg ? max_pos : max_neg;
                lc = 1 / abs_max;
                l = lc;
                nu = nu / 2;
            }
            l = nu * l;
        }

        cnt++;

        if (Sd < S)
        {
            S = Sd;

            x = xd;
            r = rd;
            J = jacobian(x, params);

            nfJ = nfJ + 1;
            A = J.transpose() * J;
            v = J.transpose() * r;
        }

		//cout << "r=" << "\n";
        //std::cout << r.norm() << "\n";
    }

    result.J = J;
    result.r = r;
    result.stopping_criteria = ERROR;
    if (cnt >= MaximumIterationNumber) result.stopping_criteria = MAXIMUM_ITERATION_REACHED;
    if (!(any_greater_than(d, epsx) && (any_greater_than(r, epsy)))) result.stopping_criteria = THRESHOLD_REACHED;
    result.no_of_iterations = cnt;
            
    return x;
}

/*
 *****************************************************
 * Simulated annealing optimizer
 *****************************************************
 */

/*VectorXd simulated_annealing(VectorXd (*fn)(VectorXd, void*), MatrixXd (*jacobian)(VectorXd, void*), void* params, VectorXd x, const double TolX, const double TolY, optimizer_result &result)
{
    const double lower_bound = -10;
    const double upper_bound = 10;
    
    VectorXd best_sol;
    double best_norm = DBL_MAX;
    optimizer_result best_result;
    
    for (size_t k = 0; k < 100; k++)
    {
        VectorXd x0 = x;
        
        for (size_t i = 0; i < x0.size(); i++)
        {
            std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
            std::default_random_engine re;
   
            x0(i) += unif(re);
        }
        
        optimizer_result local_result;
        VectorXd sol_vec = levenberg_marquardt(fn, jacobian, params, x, TolX, TolY, local_result);
        
        double r_norm = local_result.r.norm();
        if (r_norm < best_norm)
        {
            best_sol = sol_vec;
            best_norm = r_norm;
            best_result = local_result;
        }        
    }
    
    result = best_result;
    
    return best_sol;
}*/

