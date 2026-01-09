#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <fstream>

// int generateLinSystemCSR(const double L
//                         ,const int NX
//                         ,const int NY
//                         ,int const* cooRow
//                         ,int const* cooCol
//                         ,double const* cooMat
//                         ,int* csrRow
//                         // ,int* csrCol
//                         // ,double* csrMat
//                         ,double* rhs
//                         ,int& nnz 
//                         ,std::ofstream& LOG_FILE);

int coo2csr(
    const int N,
    int* cooRow,
    int* cooCol,
    double* cooMat,
    int& nnz,
    int* csrRow);

int coo2csr_sorted(int N,          // number of rows
                   int nnz,       // number of nonzeros
                   int* cooRow,   // input COO rows
                   int* cooCol,   // input COO cols (will become CSR col)
                   double* cooMat,// input COO values (will become CSR values)
                   int* csrRow);

double computeResidualCSR(int nnz,
                          int N,
                          const int* csrRow,
                          const int* csrCol,
                          const double* csrMat,
                          const double* sol,
                          const double* rhs,
                          double* const aux);


int solveJacobi2D_B(const int nnz
                   ,const int N
                   ,const int * const csrRow
                   ,const int * const csrCol
                   ,const double* const csrMat
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


// Debugging functions
void printCOO(const int nnz
              ,const int * const cooRow
              ,const int * const cooCol
              ,const double * const cooMat);
void printCSR(const int N
              ,const int * const csrRow
              ,const int * const csrCol
              ,const double * const csrMat);

#endif
