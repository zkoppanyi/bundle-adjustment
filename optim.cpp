//#define NDEBUG
#define EIGEN_NO_DEBUG 

#include <Eigen/IterativeLinearSolvers>
        
#include "optim.h"
#include "structs.h"
#include <cmath>
#include <iostream>
#include "time.h"

using namespace std;
using namespace Eigen;

bool any_greater_than(VectorXd v1, VectorXd v2)
{
    for (int ik = 0; ik < v1.size(); ik++)
    {
        if (fabs(v1(ik)) >= fabs(v2(ik)))
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

//#define PRINT_ITERATION_DETAILS

//#define SOLVER_TYPE 1  // Sparse Choelsky without analyzePatter
//#define SOLVER_TYPE 3  // Conjugate gradient; preconditioing might be needed
//#define SOLVER_TYPE 2  // Sparse Choelsky with analyzePattern 
#define SOLVER_TYPE 41 // Shur complement

typedef Eigen::Triplet<double> T;

VectorXd levenberg_marquardt(VectorXd (*fn)(VectorXd, void*), int (*jacobian)(VectorXd, Eigen::SparseMatrix<double>&, void*), void* params, VectorXd x, const double TolX, const double TolY, optimizer_result &result)
{
    const problem* prob = (problem*)params;
	const unsigned int MaximumIterationNumber = 50;

	// Residual at starting point
    VectorXd r = fn(x, params);
    double S = r.dot(r);

    size_t lr = r.size();
    size_t lx = x.size();

    double jepsx = TolX;

    Eigen::SparseMatrix<double> J(prob->sum_obs, prob->sum_unknowns);  
    jacobian(x, J, params);    

    //Debug.WriteLine("Jacobian: " + J);

    int nfJ = 2;
    Eigen::SparseMatrix<double> A = J.transpose() * J;
    VectorXd v = J.transpose() * r;

    // Automatic scaling
    MatrixXd Dd = A.diagonal().asDiagonal();   // TODO: implement to be sparse!
    size_t Ddim = min(Dd.rows(), Dd.cols());
    for (int i = 0; i < Ddim; i++)
    {
        if (Dd(i, i) == 0) Dd(i, i) = 1.0;
    }
    Eigen::SparseMatrix<double> D = Dd.sparseView();

    double Rlo = 0.25;
    double Rhi = 0.75;
    double l = 1;
    double lc = 0.75;
    unsigned int cnt = 0;

    VectorXd epsx = VectorXd::Ones(lx) * TolX;
    VectorXd epsy = VectorXd::Ones(lr) * TolY;
    VectorXd d = VectorXd::Ones(lx) * TolX;   
    
    SimplicialLLT<SparseMatrix<double>> solver; 
    bool is_pattern_anal = false;
    if (SOLVER_TYPE == 2)
    {       
        solver.analyzePattern(A);
        bool is_pattern_anal = true;
    }
    
    double res_prev = 10000;
    double dr = 1000;    
    
    while ((cnt < MaximumIterationNumber) && any_greater_than(d, epsx) && (any_greater_than(r, epsy)) && (dr > TolY))
    {       
        #ifdef PRINT_ITERATION_DETAILS
            cout << "Iteration #" << cnt << endl;
            char buffer[250];
            sprintf(buffer, "   Jacobian size: %d x %d (%d = %.1f MB)", J.rows(), J.cols(), J.rows()*J.cols(), J.rows()*J.cols()*sizeof(double)/1000000.0);
            cout << buffer << endl;
            sprintf(buffer, "   Normal matrix size: %d x %d (%d = %.1f MB)", A.rows(), A.cols(), A.rows()*A.cols(), A.rows()*A.cols()*sizeof(double)/1000000.0);
            cout << buffer << endl;
                    
            clock_t t = clock();
            clock_t t0 = t;
            clock_t tsys = t;
        #endif       
        
        Eigen::SparseMatrix<double> Av = A + l*D;               
        
        // Sparse Choelsky without analyzePatter
        if (SOLVER_TYPE == 1)
        {
            #ifdef PRINT_ITERATION_DETAILS
                cout << "   Solver type: Sparse Choelsky without pre-pattern analysis" << endl;
            #endif
                
            d = solver.compute(Av).solve(v);
        }        

        // Sparse Choelsky with analyzePattern
        if (SOLVER_TYPE == 2)
        {
            #ifdef PRINT_ITERATION_DETAILS
                cout << "   Solver type: Sparse Choelsky with pre-pattern analysis" << endl;
            #endif

            solver.factorize(Av);
            d = solver.solve(v);
        }

        // Conjugate gradient; preconditioing might be needed
        if (SOLVER_TYPE == 3)
        {
            #ifdef PRINT_ITERATION_DETAILS
                cout << "   Solver type: Conjugate gradient" << endl;
            #endif

            ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
            d = solver.compute(Av).solve(v);
        }
        
        // Schur complement       
        if (SOLVER_TYPE >= 40)
        {
            #ifdef PRINT_ITERATION_DETAILS
                cout << "   Solver type: Schur complement" << endl;
            #endif

                
            size_t n_cblock = prob->sum_unknowns - prob->start_idx_obj_pts;
            Eigen::SparseMatrix<double> Cd = Av.block(prob->start_idx_obj_pts, prob->start_idx_obj_pts, n_cblock, n_cblock); // sparse
            MatrixXd E = Av.block(0, prob->start_idx_obj_pts, prob->start_idx_obj_pts, n_cblock);
            Eigen::SparseMatrix<double> B = Av.block(0, 0, prob->start_idx_obj_pts, prob->start_idx_obj_pts); // sparse

            #ifdef PRINT_ITERATION_DETAILS
            cout << "   Initalization [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
            t = clock();
            #endif

            /*MatrixXd Cd = Av.block(prob->start_idx_obj_pts, prob->start_idx_obj_pts, n_cblock, n_cblock); // sparse
            MatrixXd Cinvd = MatrixXd::Zero(Cd.rows(), Cd.cols());
            for (size_t k = 0; k < n_cblock; k+=3)
            {
                Cinvd.block<3,3>(k, k) = Cd.block<3,3>(k, k).inverse();            
            }
            Eigen::SparseMatrix<double> Cinv = Cinvd.sparseView();*/
            
            std::vector< Eigen::Triplet<double> > tripletList;
            tripletList.reserve(n_cblock*3*3);
            for (size_t k = 0; k < n_cblock; k+=3)
            {
                Matrix3d cInv = ((Matrix3d)Cd.block(k, k, 3, 3)).inverse();                      
                tripletList.push_back(T(k+0,k+0,cInv(0,0)));
                tripletList.push_back(T(k+0,k+1,cInv(0,1)));
                tripletList.push_back(T(k+0,k+2,cInv(0,2)));
                tripletList.push_back(T(k+1,k+0,cInv(1,0)));
                tripletList.push_back(T(k+1,k+1,cInv(1,1)));
                tripletList.push_back(T(k+1,k+2,cInv(1,2)));
                tripletList.push_back(T(k+2,k+0,cInv(2,0)));
                tripletList.push_back(T(k+2,k+1,cInv(2,1)));
                tripletList.push_back(T(k+2,k+2,cInv(2,2)));                
            }
            Eigen::SparseMatrix<double> Cinv( Cd.rows(), Cd.cols() );
            Cinv.setFromTriplets(tripletList.begin(), tripletList.end());

            #ifdef PRINT_ITERATION_DETAILS
                cout << "   C inverse[s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
                t = clock();
            #endif

            VectorXd v2 = v.head(prob->start_idx_obj_pts);
            VectorXd w = v.tail(n_cblock);        
            
            SparseMatrix<double> Sch = B - E*Cinv*E.transpose(); // this seems to be very computation intensive here
            VectorXd rhs = v2 - E*Cinv*w;                

            // Solve the reduced camera system

            // dense version
            //VectorXd dy =(B - E*Cinv*E.transpose()).llt().solve(v2 - E*Cinv*w); 

            // sparse version            
            if (!is_pattern_anal)
            {
                 solver.analyzePattern(Sch);
                 is_pattern_anal = true;
            }
            solver.factorize(Sch);
            VectorXd dy = solver.solve(rhs);
            
            //VectorXd dy = solver.compute(Sch).solve(rhs);

            VectorXd dz = Cinv*(w-E.transpose()*dy);        
            d << dy, dz;
        }
               

        #ifdef PRINT_ITERATION_DETAILS
        cout << "   Chloesky [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif        
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   dz calculation [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        cout << "   Normal system time [s]: " << (float)(clock() - tsys)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif
        
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
            {
                nu = 2;
            }
            else if (nu > 10)
            {
                nu = 10;
            }

            if (l == 0)
            {
                VectorXd diag =  MatrixXd(A).inverse().diagonal();
                double max_pos = diag.maxCoeff();
                double max_neg = abs(diag.minCoeff());
                
                /*VectorXd diag = A.diagonal();
                double max_pos = diag(0);
                double max_neg = diag(0);
                for (int k=0; k<diag.size();k++)
                {
                    double diagk = diag(k) == 0 ? 0 : diag(k);
                	max_pos = max_pos > 1/diagk ? max_pos : 1/diagk;
                	max_neg = max_neg < 1/diagk ? max_neg : 1/diagk;
                }
				max_neg = abs(max_neg);*/

                double abs_max = max_pos > max_neg ? max_pos : max_neg;
                lc = 1 / abs_max;
                l = lc;
                nu = nu / 2;
            }
            l = nu * l;
        }

        cnt++;
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   Adjusting time [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif

        if (Sd < S)
        {
            S = Sd;

            x = xd;
            r = rd;
            jacobian(x, J, params);
            
            #ifdef PRINT_ITERATION_DETAILS
            cout << "   Populate Jacobian   [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
            cout << "   # of non-zeros in J [s]: " << J.nonZeros() << endl;
            t = clock();
            #endif

            nfJ = nfJ + 1;
            A = J.transpose() * J;
            v = J.transpose() * r;
            
            #ifdef PRINT_ITERATION_DETAILS
            cout << "   Update A and v [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
            t = clock();
            #endif
        }
       
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   Rest of iteration [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
		//cout << "r=" << "\n";
        std::cout << "   Residual : " << r.norm() << "\n";
        cout << "   Total iteration [s]: " << (float)(clock() - t0)/CLOCKS_PER_SEC << endl;
        //std::cout << std::endl;
        //std::cout << d << std::endl;
        #endif
        
        double r_norm = r.norm();
        dr =  abs(r_norm - res_prev);
        res_prev = r_norm;
        
        std::cout << "#" << cnt << "   Residual : " << r_norm << " [" << dr << "]"  "\n";
        //printf("-- LOG: #i Residual: %.6f [%.6f]\n", cnt, r_norm, dr);
    }

    result.J = J;
    result.r = r;
    result.stopping_criteria = ERROR;
    if (cnt >= MaximumIterationNumber) result.stopping_criteria = MAXIMUM_ITERATION_REACHED;
    if (!(any_greater_than(d, epsx) && (any_greater_than(r, epsy)) && (dr > TolY))) result.stopping_criteria = THRESHOLD_REACHED;
    result.no_of_iterations = cnt;
            
    return x;
}


/*VectorXd levenberg_marquardt(VectorXd (*fn)(VectorXd, void*), int (*jacobian)(VectorXd, Eigen::SparseMatrix<double>&, void*), void* params, VectorXd x, const double TolX, const double TolY, optimizer_result &result)
{
    const problem* prob = (problem*)params;
	const unsigned int MaximumIterationNumber = 100;

	// Residual at starting point
    VectorXd r = fn(x, params);
    double S = r.dot(r);

    size_t lr = r.size();
    size_t lx = x.size();

    double jepsx = TolX;

    Eigen::SparseMatrix<double> J(prob->sum_obs, prob->sum_unknowns);    
    jacobian(x, J, params);
    //Debug.WriteLine("Jacobian: " + J);

    int nfJ = 2;
    MatrixXd A = MatrixXd(J.transpose() * J);
    VectorXd v = MatrixXd(J.transpose() * r);

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
        d = (A + l*D).ldlt().solve(v); // use ldlt for stability.

        
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

            jacobian(x, J, params);


            nfJ = nfJ + 1;
            A = MatrixXd(J.transpose() * J);
            v = MatrixXd(J.transpose() * r);
        }

        //std::cout << cnt << ".   r= " << r.norm() << " d = " << d.norm() << "\n";
        printf("%i. r = %.6f d = %.6f\n", cnt, r.norm(), d.norm());
    }

    result.J = J;
    result.r = r;
    result.stopping_criteria = ERROR;
    if (cnt >= MaximumIterationNumber) result.stopping_criteria = MAXIMUM_ITERATION_REACHED;
    if (!(any_greater_than(d, epsx) && (any_greater_than(r, epsy)))) result.stopping_criteria = THRESHOLD_REACHED;
    result.no_of_iterations = cnt;
            
    return x;
}*/

/*VectorXd levenberg_marquardt(VectorXd (*fn)(VectorXd, void*), int (*jacobian)(VectorXd, Eigen::SparseMatrix<double>&, void*), void* params, VectorXd x, const double TolX, const double TolY, optimizer_result &result)
{
    const problem* prob = (problem*)params;
	const unsigned int MaximumIterationNumber = 50;

	// Residual at starting point
    VectorXd r = fn(x, params);
    double S = r.dot(r);

    size_t lr = r.size();
    size_t lx = x.size();

    double jepsx = TolX;

    Eigen::SparseMatrix<double> J(prob->sum_obs, prob->sum_unknowns);    

    jacobian(x, J, params);    

    //Debug.WriteLine("Jacobian: " + J);

    int nfJ = 2;
    Eigen::SparseMatrix<double> A = J.transpose() * J;
    VectorXd v = J.transpose() * r;

    // Automatic scaling
    MatrixXd Dd = A.diagonal().asDiagonal();   // TODO: implement to be sparse!
    size_t Ddim = min(Dd.rows(), Dd.cols());
    for (int i = 0; i < Ddim; i++)
    {
        if (Dd(i, i) == 0) Dd(i, i) = 1.0;
    }
    Eigen::SparseMatrix<double> D = Dd.sparseView();

    double Rlo = 0.25;
    double Rhi = 0.75;
    double l = 1;
    double lc = 0.75;
    unsigned int cnt = 0;

    VectorXd epsx = VectorXd::Ones(lx) * TolX;
    VectorXd epsy = VectorXd::Ones(lr) * TolY;

    VectorXd d = VectorXd::Ones(lx) * TolX;

    SimplicialLLT<SparseMatrix<double>> solver; 
    solver.analyzePattern(A);
    
    //while (cnt < MaximumIterationNumber)
    while ((cnt < MaximumIterationNumber) && any_greater_than(d, epsx) && (any_greater_than(r, epsy)))
    {
        // negative solution increment
        //d = SolveLinearEquationSystem(A.Add((l.Multiply(D))), v);
        //d = MatrixXd((A + l*D)).ldlt().solve(v); // use ldlt for stability.
        //d = (A + l*D).inverse()*(v); 
        //d = (A + l*D).llt().solve(v); 
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "Iteration #" << cnt << endl;
        clock_t t = clock();
        clock_t t0 = t;
        clock_t tsys = t;
        #endif       
        
        Eigen::SparseMatrix<double> Av = A + l*D;        

        // Schur complement
        // Sparse Choelsky seems faster than all matrix manipluation under here
        /*size_t n_cblock = prob->sum_unknowns - prob->start_idx_obj_pts;
        MatrixXd Cd = Av.block(prob->start_idx_obj_pts, prob->start_idx_obj_pts, n_cblock, n_cblock); // sparse
        MatrixXd E = Av.block(0, prob->start_idx_obj_pts, prob->start_idx_obj_pts, n_cblock);
        Eigen::SparseMatrix<double> B = Av.block(0, 0, prob->start_idx_obj_pts, prob->start_idx_obj_pts); // sparse
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   Initalization [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif
        
        MatrixXd Cinvd = MatrixXd::Zero(Cd.rows(), Cd.cols());
        for (size_t k = 0; k < n_cblock; k+=3)
        {
            Cinvd.block<3,3>(k, k) = Cd.block<3,3>(k, k).inverse();            
        }
        Eigen::SparseMatrix<double> Cinv = Cinvd.sparseView();

        #ifdef PRINT_ITERATION_DETAILS
        cout << "   C inverse[s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif
        
        VectorXd v2 = v.head(prob->start_idx_obj_pts);
        VectorXd w = v.tail(n_cblock);        
        
        SparseMatrix<double> Sch = B - E*Cinv*E.transpose(); // this seems to be very computation intensive here
        VectorXd rhs = v2 - E*Cinv*w;                
        //VectorXd dy =(B - E*Cinv*E.transpose()).llt().solve(v2 - E*Cinv*w); // dense version
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   S calculation [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif
        
        SimplicialLLT<SparseMatrix<double>> solver; // sparse version
        VectorXd dy = solver.compute(Sch).solve(rhs);
         
        VectorXd dz = Cinv*(w-E.transpose()*dy);        
        d << dy, dz;*/
        
        // Sparse Choelsky without analyzePatter
        /*VectorXd d = solver.compute(Av).solve(v);
        cout << "d= " << d.norm() << endl;

        // Sparse Choelsky with analyzePattern
        //solver.factorize(Av);
        //VectorXd d = solver.solve(v);

        // Conjugate gradient; preconditioing might be needed
        //ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
        //VectorXd d = solver.compute(Av).solve(v);

        cout << "   l: " << l  << endl;
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   Chloesky [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif       
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   dz calculation [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        cout << "   Normal system time [s]: " << (float)(clock() - tsys)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif
        
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
            {
                nu = 2;
            }
            else if (nu > 10)
            {
                nu = 10;
            }

            //if (l < lc)
            if (l == 0)
            {
                VectorXd diag = MatrixXd(A).inverse().diagonal();
                //VectorXd diag = A.diagonal();
                double max_pos = diag.maxCoeff();
                double max_neg = abs(diag.minCoeff());
                /*double max_pos = diag(0);
                double max_neg = diag(0);
                for (int k=0; k<diag.size();k++)
                {
                    //double diagk = diag(k) == 0 ? 0 : 1/diag(k);
                    double diagk = diag(k);
                	max_pos = max_pos > diagk ? max_pos : diagk;
                	max_neg = max_neg < diagk ? max_neg : diagk;
                }
				max_neg = abs(max_neg);*//*

                double abs_max = max_pos > max_neg ? max_pos : max_neg;
                lc = 1 / abs_max;
                l = lc;
                nu = nu / 2;
            }
            l = nu * l;
        }

        cnt++;
        
        #ifdef PRINT_ITERATION_DETAILS
        cout << "   Adjusting time [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
        t = clock();
        #endif

        //if (Sd < S)
        {
            S = Sd;
            x = xd;
            r = rd;
            jacobian(x, J, params);
            
            #ifdef PRINT_ITERATION_DETAILS
            cout << "   Populate Jacobian   [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
            cout << "   # of non-zeros in J [s]: " << J.nonZeros() << endl;
            t = clock();
            #endif

            nfJ = nfJ + 1;
            A = J.transpose() * J;
            v = J.transpose() * r;
            
            #ifdef PRINT_ITERATION_DETAILS
            cout << "   Update A and v [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
            t = clock();
            #endif
        }
        
        std::cout << "   Residual : " << r.norm() << "\n";

        #ifdef PRINT_ITERATION_DETAILS
        cout << "   Rest of iteration [s]: " << (float)(clock() - t)/CLOCKS_PER_SEC << endl;
		//cout << "r=" << "\n";
        std::cout << "   Residual : " << r.norm() << "\n";
        cout << "   Total iteration [s]: " << (float)(clock() - t0)/CLOCKS_PER_SEC << endl;
        //std::cout << std::endl;
        //std::cout << d << std::endl;
        #endif
    }

    result.J = J;
    result.r = r;
    result.stopping_criteria = ERROR;
    if (cnt >= MaximumIterationNumber) result.stopping_criteria = MAXIMUM_ITERATION_REACHED;
    if (!(any_greater_than(d, epsx) && (any_greater_than(r, epsy)))) result.stopping_criteria = THRESHOLD_REACHED;
    result.no_of_iterations = cnt;
            
    return x;
}*/

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

