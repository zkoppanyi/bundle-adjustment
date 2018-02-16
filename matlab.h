#ifndef MATLAB_H
#define MATLAB_H

#include <Eigen/Dense>
#include "mex.h"
#include "structs.h"
#include "core.h"

#define LOG(S) mexPrintf("-- LOG: %s\n", S);
#define LOG2I(S,I) mexPrintf("-- LOG: %s %i\n", S, I);
#define LOG2F(S,F) mexPrintf("-- LOG: %s %.5f\n", S, F);

#define LOG_ERROR(S) mexPrintf("!! ERR: %s\n", S); mexErrMsgTxt("Error occured! Exiting...");
#define LOG_WARNING(S) mexPrintf("-! WRN: %s\n", S);


/*#define ASSERT( isOK )        if (!(isOK)) {  \
                                    mexPrintf("Assert failed: %s:%s (%s)\n", __FILE__, __LINE__, __FUNCTION__);  \
                                    mexErrMsgTxt("Assertation failed! Exiting...");  \
                                };*/
                                
#define ASSERT( isOK )        if (!(isOK)) {  \
                                    mexPrintf("Assert failed: %s:%s (%s)\n", __FILE__, __LINE__, __FUNCTION__);  \
                                };                                
 
#define GET(A,i,j,n) A[(i) + (j)*n]

void print(Eigen::MatrixXd A);
void print(Eigen::VectorXd v);

void eigen2mat(const Eigen::MatrixXd &A, mxArray* &ret);
void eigen2mat(double d, mxArray* &ret);

void create_stoch_struct(const optimizer_result &optimizer_result, const stochastic_params &stoch, mxArray* &ret);
void create_problem_struct(const problem &prob, mxArray* &ret);

int extract_problem_from_arguments(int nrhs, const mxArray *prhs[], problem &prob, int mode = 1 );

#endif