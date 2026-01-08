#include "functions.h"

using namespace std;

inline int ij2l(const int i, const int j, const int Nx) {
   return (j-1)*(Nx-1)+(i-1);
}

int generateLinearSystemCSR(const double L
                           ,const int NX
                           ,const int NY
                           ,int const* cooRow
                           ,int const* cooCol
                           ,double const* cooMat
                           ,int* csrRow
                           ,int* csrCol
                           ,double* csrMat
                           ,double* rhs
                           ,int& nnz 
                           ,ofstream& LOG_FILE){
   

   // REVISIT TO CHECK IF EVERYTHING IS CORRECT!!!!!!!!!!!!!!!!!!!!!!!!!!!
   for (int k=0; k < nnz; ++k){
      csrCol[k] = cooCol[k];
      csrMat[k] = cooMat[k];
   }
   // TO BE DONE: fill csrRow array
   for (int i=0; i<= (NX-1)*(NY-1); ++i){
      csrRow[i] = 0;
   }

   // count entries per row:
   for (int k = 0; k < nnz; ++k){
      csrRow[cooRow[k] + 1]++; // note the +1 shift as first entry is zero
   }

   // Prefix sum to get row pointers:
   for (int i = 0; i < (NX-1)*(NY-1); ++i){
      csrRow[i + 1] += csrRow[i];
   }

   return EXIT_SUCCESS;
}

double computeResidualCSR(int nnz,
                          int N,
                          const int* csrRow,
                          const int* csrCol,
                          const double* csrMat,
                          const double* sol,
                          const double* rhs,
                          const double* aux) {
   
   // TO BE DONE: implement residual computation for CSR format!!!!!!!!!!!
   return 0.0;
   }

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
                   ,ofstream& LOG_FILE) {

   res = 2*TOL; // initialize residual

   iters = 0; // initialize iteration counter

   while (res > TOL && iters < MAX_ITERS) {

      for (int i = 0; i<N; ++i) { // loop over all states/rows
         double sum = 0.0;
         double aii = 0.0; // diagonal entry

         for (int k = csrRow[i]; k < csrRow[i + 1]; ++k) {
            int j = csrCol[k]; // extract current column
            int val = csrMat[k]; // extract current value

            if (i == j)
               aii = val;
            else
               sum += val * sol[j]; // accumulate sum with current solution times matrix entry
         }
         aux[i] = (rhs[i] - sum) / aii; // update aux with new solution.
      }

      // TODO: COMPUTE RESIDUALCSR!!!!!!!!!!!!!!!!!!
      // compute residual, swtiching sol and aux as needed
      res = computeResidualCSR(nnz,N
                              ,csrRow,csrCol,csrMat
                              ,aux,rhs,sol); // note switch between aux and sol

      for (int l = 0; l<N; ++l) {
         sol[l] = aux[l]; // update solution with aux values
      }

      ++iters;
      if (iters == MAX_ITERS) break; // break if max iters reached
      if (iters < MAX_ITERS) {
         LOG_FILE << "Successfull convergence" << endl;
         LOG_FILE << " - residual achieved " << res << endl;
         LOG_FILE << " - iterations " << iters << endl;
         return EXIT_SUCCESS;
      } else { 
         LOG_FILE << "Failure of convergence procedure" << endl;
         LOG_FILE << " - residual achieved " << res << endl;
         LOG_FILE << " - iterations " << iters << endl;
         return EXIT_FAILURE;
   }
   }
   

   //return EXIT_SUCCESS;
}

int solveJacobi2D_C(const double L
                   ,const int NX, const int NY
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,const double* const rhs
                   ,double& res , int& iters
                   ,double* const aux
                   ,std::ofstream& LOG_FILE) {

   return EXIT_SUCCESS;
}

