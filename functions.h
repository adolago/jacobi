#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <fstream>

int solveJacobi2D_B(const int nnz
                   ,const int N
                   ,const int * const cooRow
                   ,      int * const csrRow
                   ,const int * const cooCol
                   ,const double * const cooMat
                   ,const double* const rhs
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,double& res , int& iters
                   ,double* const aux
                   ,std::ofstream& LOG_FILE);

int solveJacobi2D_C(const double L
                   ,const int NX, const int NY
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,const double* const rhs
                   ,double& res , int& iters
                   ,double* const aux
                   ,std::ofstream& LOG_FILE);

#endif
