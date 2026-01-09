#include "functions.h"
#include <iostream>

using namespace std;

inline int ij2l(const int i, const int j, const int Nx) {
   return (j-1)*(Nx-1)+(i-1);
}

#include <algorithm>
#include <vector>
#include <numeric>

int coo2csr(
    const int N,
    int* cooRow,
    int* cooCol,
    double* cooMat,
    int& nnz,
    int* csrRow) {

   // sort COO entries by (row, col)
   std::vector<int> perm(nnz);
   std::iota(perm.begin(), perm.end(), 0);

   std::sort(perm.begin(), perm.end(),
      [&](int a, int b) {
         if (cooRow[a] != cooRow[b])
               return cooRow[a] < cooRow[b];
         return cooCol[a] < cooCol[b];
      });

   // Apply permutation in-place using temporaries
   std::vector<int> rowTmp(nnz), colTmp(nnz);
   std::vector<double> valTmp(nnz);

   for (int k = 0; k < nnz; ++k) {
      rowTmp[k] = cooRow[perm[k]];
      colTmp[k] = cooCol[perm[k]];
      valTmp[k] = cooMat[perm[k]];
   }

   std::copy(rowTmp.begin(), rowTmp.end(), cooRow);
   std::copy(colTmp.begin(), colTmp.end(), cooCol);
   std::copy(valTmp.begin(), valTmp.end(), cooMat);

   // Compress duplicates (sum stencil contributions)
   int write = 0;

   for (int read = 0; read < nnz; ++read) {
      if (write == 0 ||
         cooRow[read] != cooRow[write - 1] ||
         cooCol[read] != cooCol[write - 1]) {

         cooRow[write] = cooRow[read];
         cooCol[write] = cooCol[read];
         cooMat[write] = cooMat[read];
         ++write;
      } else {
         // same (row,col) → sum contributions
         cooMat[write - 1] += cooMat[read];
      }
   }

   nnz = write;  // compressed nnz

   // Build CSR row pointers

   std::fill(csrRow, csrRow + N + 1, 0);

   for (int k = 0; k < nnz; ++k) {
      csrRow[cooRow[k] + 1]++;
   }

   for (int i = 0; i < N; ++i) {
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
                          double* const aux) {
   
   // set aux to zero
   for (int l = 0; l<N; ++l) {
      aux[l] = 0.0;
   }

   // aux = A*sol
   for (int i = 0; i<N; ++i) {
      for (int k = csrRow[i]; k < csrRow[i + 1]; ++k) {
         aux[i] += csrMat[k] * sol[csrCol[k]]; // accumulate row i by multiplying row i of A with sol
      }
   }
   // compute residual
   double res = 0.0;
   for (int i = 0; i<N; ++i) {
      res += sqrt(pow(rhs[i] - aux[i], 2));
   }
   res /= sqrt(N);
   return res;
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
            double val = csrMat[k]; // extract current value

            if (i == j){
               aii = val;
               // std::cout << "Diagonal entry at row " << i << " : " << aii << std::endl;
               if (aii == 0.0){
                  // std::cout << "Zero diagonal entry found at row " << i << ". Exiting." << std::endl;
                  return EXIT_FAILURE;
               }
            }
            else
               sum += val * sol[j]; // accumulate sum with current solution times matrix entry
            
         }
         aux[i] = (rhs[i] - sum) / aii; // update aux with new solution.
      }

      // compute residual, swtiching sol and aux as needed
      res = computeResidualCSR(nnz,N
                              ,csrRow,csrCol,csrMat
                              ,aux,rhs,sol); // note switch between aux and sol
      // std::cout << "Residual at iteration " << iters << " : " << res << std::endl;
      for (int l = 0; l<N; ++l) {
         sol[l] = aux[l]; // update solution with aux values
      }

      ++iters;
      if (iters == MAX_ITERS) break; // break if max iters reached
   }

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

void printCOO(const int nnz
              ,const int * const cooRow
              ,const int * const cooCol
              ,const double * const cooMat) {
   cout << "COO matrix format: " << endl;
   for (int k = 0; k<nnz; ++k) {
      cout << "(" << cooRow[k] << "," << cooCol[k] << ") : " << cooMat[k] << endl;
   }
}
void printCSR(const int N
              ,const int * const csrRow
              ,const int * const csrCol
              ,const double * const csrMat) {
   cout << "CSR matrix format: " << endl;
   for (int i = 0; i<N; ++i) {
      cout << "Row " << i << " : ";
      for (int k = csrRow[i]; k<csrRow[i+1]; ++k) {
         cout << "(" << csrCol[k] << ") : " << csrMat[k] << " ; ";
      }
      cout << endl;
   }
}