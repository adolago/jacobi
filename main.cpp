#include <cstdlib>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

#if defined(_WIN32) || defined(_WIN64)
const string PATH_SEP = "\\";
#else
const string PATH_SEP = "/";
#endif

inline int ij2l(const int i, const int j, const int Nx) {
   return (j-1)*(Nx-1)+(i-1);
}

int generateLinSystemCOO(const double L
                        ,const int NX
                        ,const int NY
                        ,int* const cooRow
                        ,int* const cooCol
                        ,double* const cooMat
                        ,double* const rhs
                        ,int& nnz 
                        ,ofstream& LOG_FILE);

double computeError(const double L
                   ,const int NX
                   ,const int NY
                   ,const double* const sol
                   ,ofstream& LOG_FILE);

int solveJacobi2D_A(const int nnz
                   ,const int N
                   ,const int * const cooRow
                   ,const int * const cooCol
                   ,const double * const cooMat
                   ,const double* const rhs
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,double& res , int& iters
                   ,double* const aux
                   ,ofstream& LOG_FILE);

double computeResidual(const int nnz
                   ,const int N
                   ,const int * const cooRow
                   ,const int * const cooCol
                   ,const double * const cooMat
                   ,const double * const sol 
                   ,const double * const rhs 
                   ,      double * const aux);

#include "functions.h"

int main(int argc, char* argv[]) {
  
   constexpr static const int SIZE_STENCIL = 5;
   constexpr static const double L = 1.0;

   const int NX = atoi(argv[1]);
   const double TOL = atof(argv[2]);
   const int MAX_ITERS = atoi(argv[3]);
   const string OUTPUT_DIR = argv[4];

   ofstream LOG_FILE(  OUTPUT_DIR + PATH_SEP + "screen.log");

   const string fnmA = OUTPUT_DIR + PATH_SEP + "results_JacobiA.bin";
   FILE* fidA = fopen(fnmA.c_str(),"wb");
   const string fnmB = OUTPUT_DIR + PATH_SEP + "results_JacobiB.bin";
   FILE* fidB = fopen(fnmB.c_str(),"wb");
   const string fnmC = OUTPUT_DIR + PATH_SEP + "results_JacobiC.bin";
   FILE* fidC = fopen(fnmC.c_str(),"wb");

   const int NY = NX;

   const int N = (NX-1)*(NY-1);

   const int NNZ_MAX = N*SIZE_STENCIL;


   int ierr = EXIT_SUCCESS;

   int*    const cooRow = new    int[NNZ_MAX];
   int*    const csrRow = new    int[NNZ_MAX];
   int*    const cooCol = new    int[NNZ_MAX];
   double* const cooMat = new double[NNZ_MAX];
   double* const rhs = new double[N];
   double* const sol = new double[N];
   double* const aux = new double[N];
   int nnz; // to be defined after linear system is built

   ierr = ierr | generateLinSystemCOO(L,NX,NY
                                     ,cooRow,cooCol,cooMat
                                     ,rhs
                                     ,nnz
                                     ,LOG_FILE);

   int iters;
   double res;
   double err;
   double seconds; 
   chrono::steady_clock::time_point tp0;
   chrono::duration<double> tDuration;

   // solver Jacobi A --------------------------------------- //
   LOG_FILE << endl << "SOLVER JACOBI A" << endl << endl;

   // reset solution
   for (int i=0;i<N;++i) sol[i] = 0.5;

   // solve
   tp0 = chrono::steady_clock::now();
   ierr = ierr | solveJacobi2D_A(nnz,N
                                ,cooRow,cooCol,cooMat
                                ,rhs
                                ,TOL,MAX_ITERS
                                ,sol,res,iters
                                ,aux
                                ,LOG_FILE);
   tDuration = chrono::steady_clock::now()-tp0;

   err  = computeError(L,NX,NY,sol,LOG_FILE);

   seconds = tDuration.count();

   // print results
   LOG_FILE << "iters   : " << iters << endl;
   LOG_FILE << "res     : " << res   << endl;
   LOG_FILE << "err     : " << err   << endl;
   LOG_FILE << "seconds : " << seconds << endl;
   
   // store results results
   fwrite(&iters   , sizeof(int)   , 1, fidA);
   fwrite(&res     , sizeof(double), 1, fidA);
   fwrite(&err     , sizeof(double), 1, fidA);
   fwrite(&seconds , sizeof(double), 1, fidA);

   // solver Jacobi B --------------------------------------- //
   LOG_FILE << endl << "SOLVER JACOBI B" << endl << endl;

   // reset solution
   for (int i=0;i<N;++i) sol[i] = 0.5;

   // solve
   tp0 = chrono::steady_clock::now();
   ierr = ierr | solveJacobi2D_B(nnz,N
                                ,cooRow,csrRow,cooCol,cooMat
                                ,rhs
                                ,TOL,MAX_ITERS
                                ,sol,res,iters
                                ,aux
                                ,LOG_FILE);
   tDuration = chrono::steady_clock::now()-tp0;

   err  = computeError(L,NX,NY,sol,LOG_FILE);

   seconds = tDuration.count();

   // print results
   LOG_FILE << "iters : " << iters << endl;
   LOG_FILE << "res   : " << res   << endl;
   LOG_FILE << "err   : " << err   << endl;
   LOG_FILE << "wall t: " << seconds << endl;
   
   // store results results
   fwrite(&iters   , sizeof(int)   , 1, fidB);
   fwrite(&res     , sizeof(double), 1, fidB);
   fwrite(&err     , sizeof(double), 1, fidB);
   fwrite(&seconds , sizeof(double), 1, fidB);
 
   // solver Jacobi C --------------------------------------- //
   LOG_FILE << endl << "SOLVER JACOBI C" << endl << endl;

   // reset solution
   for (int i=0;i<N;++i) sol[i] = 0.5;

   // solve
   tp0 = chrono::steady_clock::now();
   ierr = ierr | solveJacobi2D_C(L,NX,NY
                                ,TOL,MAX_ITERS
                                ,sol,rhs,res,iters
                                ,aux
                                ,LOG_FILE);
   tDuration = chrono::steady_clock::now()-tp0;

   err  = computeError(L,NX,NY,sol,LOG_FILE);

   seconds = tDuration.count();

   // print results
   LOG_FILE << "iters : " << iters << endl;
   LOG_FILE << "res   : " << res   << endl;
   LOG_FILE << "err   : " << err   << endl;
   LOG_FILE << "wall t: " << seconds << endl;
   
   // store results results
   fwrite(&iters   , sizeof(int)   , 1, fidB);
   fwrite(&res     , sizeof(double), 1, fidB);
   fwrite(&err     , sizeof(double), 1, fidB);
   fwrite(&seconds , sizeof(double), 1, fidB);
   LOG_FILE.close();

   fclose(fidA);
   fclose(fidB);
   fclose(fidC);

   delete[] cooRow;
   delete[] cooCol;
   delete[] cooMat;
   delete[] rhs;
   delete[] sol;
   delete[] aux;

   return EXIT_SUCCESS;

};

int generateLinSystemCOO(const double L
                        ,const int NX
                        ,const int NY
                        ,int* const cooRow
                        ,int* const cooCol
                        ,double* const cooMat
                        ,double* const rhs
                        ,int& nnz
                        ,ofstream& LOG_FILE){

   const double DX  = L/static_cast<double>(NX);
   const double DY  = L/static_cast<double>(NY);
   const double DXSQ  = pow(DX,2);
   const double DYSQ  = pow(DY,2);


   const double a10 = 1.0/DYSQ;
   const double a01 = 1.0/DXSQ;
   const double a11 = -(2.0/DXSQ+2.0/DYSQ);
   const double a21 = 1.0/DXSQ;
   const double a12 = 1.0/DYSQ;
  
   // fill matrix
   //LOG_FILE << "main nested loop " << endl;
   nnz=0;
   for (int j=2; j<NY-1; ++j) {
      for (int i=2; i<NX-1; ++i) {

         const int l = ij2l(i,j,NX);

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;
         cooRow[nnz+3] = l;
         cooRow[nnz+4] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l-(NX-1);
         cooCol[nnz+2] = l-1;
         cooCol[nnz+3] = l+1;
         cooCol[nnz+4] = l+(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a10;
         cooMat[nnz+2] = a01;
         cooMat[nnz+3] = a21;
         cooMat[nnz+4] = a12;

         nnz+=5;

         rhs[l] = 0.0;

         
      }
   }

   // Exceptions with west boundary
   //LOG_FILE << "West boundary " << endl;
   {const int i = 1;
      for (int j=2; j<NY-1; ++j) {

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;
         cooRow[nnz+3] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l-(NX-1);
         cooCol[nnz+2] = l+1;
         cooCol[nnz+3] = l+(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a10;
         cooMat[nnz+2] = a21;
         cooMat[nnz+3] = a12;

         nnz+=4;

         rhs[l] = 0.0;
         
      }
   }


   // Exceptions with east boundary
   //LOG_FILE << "East boundary " << endl;
   { const int i = NX-1;
      for (int j=2; j<NY-1; ++j) {

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;
         cooRow[nnz+3] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l-(NX-1);
         cooCol[nnz+2] = l-1;
         cooCol[nnz+3] = l+(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a10;
         cooMat[nnz+2] = a01;
         cooMat[nnz+3] = a12;

         nnz+=4;

         rhs[l] = 0.0;
         
      }
   }

   // Exceptions with south boundary
   //LOG_FILE << "South boundary " << endl;
   {const int j = 1;
      for (int i=2; i<NX-1; ++i) {

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;
         cooRow[nnz+3] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l-1;
         cooCol[nnz+2] = l+1;
         cooCol[nnz+3] = l+(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a01;
         cooMat[nnz+2] = a21;
         cooMat[nnz+3] = a12;

         nnz+=4;

         rhs[l] = 0.0;
         
      }
   }

   // Exceptions with north boundary
   //LOG_FILE << "North boundary " << endl;
   //LOG_FILE << "a12 = " << a12 << endl;
   {const int j = NY-1;
      for (int i=2; i<NX-1; ++i) {

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;
         cooRow[nnz+3] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l-(NX-1);
         cooCol[nnz+2] = l-1;
         cooCol[nnz+3] = l+1;

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a10;
         cooMat[nnz+2] = a01;
         cooMat[nnz+3] = a21;

         nnz+=4;

         const double xi = static_cast<double>(i)*DX;
         //LOG_FILE << "xi = " << xi << endl;
         rhs[l] = -a12*std::sin(M_PI*xi);
         //LOG_FILE << "rhs = " << rhs[l] << endl;
         
      }
   }

   // Exceptions with southwest boundary
   //LOG_FILE << "Southwest boundary " << endl;
   {const int i = 1;
      {const int j = 1;

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l+1;
         cooCol[nnz+2] = l+(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a21;
         cooMat[nnz+2] = a12;

         nnz+=3;

         rhs[l] = 0.0;
         
      }
   }

   // Exceptions with southeast boundary
   //LOG_FILE << "SouthEast boundary " << endl;
   {const int i = NX-1;
      {const int j = 1;

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l-1;
         cooCol[nnz+2] = l+(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a01;
         cooMat[nnz+2] = a12;

         nnz+=3;

         rhs[l] = 0.0;
         
      }
   }

   // Exceptions with northwest boundary
   //LOG_FILE << "Northwest boundary " << endl;
   {const int i = 1;
      {const int j = NY-1;

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l+1;
         cooCol[nnz+2] = l-(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a21;
         cooMat[nnz+2] = a10;

         nnz+=3;

         const double xi = static_cast<double>(i)*DX;
         rhs[l] = -a12*std::sin(M_PI*xi);
         
      }
   }

   // Exceptions with northeast boundary
   //LOG_FILE << "Northeast boundary " << endl;
   {const int i = NX-1;
      {const int j = NY-1;

         const int l = ij2l(i,j,NX);

         //LOG_FILE << "("<<i<<","<<j<<") --> " << l << endl;

         cooRow[nnz+0] = l;
         cooRow[nnz+1] = l;
         cooRow[nnz+2] = l;

         cooCol[nnz+0] = l;
         cooCol[nnz+1] = l-1;
         cooCol[nnz+2] = l-(NX-1);

         cooMat[nnz+0] = a11;
         cooMat[nnz+1] = a01;
         cooMat[nnz+2] = a10;

         nnz+=3;

         const double xi = static_cast<double>(i)*DX;
         rhs[l] = -a12*std::sin(M_PI*xi);
         
      }
   }

   return EXIT_SUCCESS;
}


int solveJacobi2D_A(const int nnz
                   ,const int N
                   ,const int * const cooRow
                   ,const int * const cooCol
                   ,const double* const cooMat
                   ,const double* const rhs
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,double& res , int& iters
                   ,double* const aux
                   ,ofstream& LOG_FILE) {

   res = 2*TOL;

   iters = 0;
   while (res>TOL) {
      int inz = 0;
      while (inz <= nnz-1 ) {
         double sum = 0.0;
         double a11;
         const int l = cooRow[inz];
         int row=cooRow[inz];
         while ( row == l and inz <= nnz-1) {
            const int col=cooCol[inz];
            const double coeff=cooMat[inz];
            if (row == col) { 
               a11 = coeff;
            } else {
               sum+= coeff*sol[col];
            }
            inz+=1;
            row=cooRow[inz];
         }
         aux[l] = (rhs[l]-sum)/a11;
      }

      // compute residual
      res = computeResidual(nnz,N
                           ,cooRow,cooCol,cooMat
                           ,aux,rhs,sol);

      for (int l = 0; l<N;++l) {
         sol[l] = aux[l];
      }

      ++iters;
      if (iters == MAX_ITERS) break;
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


double computeResidual(const int nnz
                   ,const int N
                   ,const int * const cooRow
                   ,const int * const cooCol
                   ,const double * const cooMat
                   ,const double * const sol 
                   ,const double * const rhs 
                   ,      double * const aux) {

   for (int l = 0; l<N;++l) {
      aux[l] = 0.0;
   }

   for (int l = 0; l<nnz;++l) {
      const int row=cooRow[l];
      const int col=cooCol[l];
      const double coeff = cooMat[l];
      aux[row]+=coeff*sol[col];

   }

   double res = 0.0;
   for (int l = 0; l<N;++l) {
      res += sqrt(pow(rhs[l]-aux[l],2));
      //res = max(res,abs(rhs[l]-aux[l]));
   }
   res/=sqrt(N);
   return res;
}

double computeError(const double L
                   ,const int NX
                   ,const int NY
                   ,const double* const sol
                   ,ofstream& LOG_FILE){

   const double DX  = L/static_cast<double>(NX);
   const double DY  = L/static_cast<double>(NY);


   double err = 0.0;
   for (int j=1; j<=NY-1; ++j) {
      const double yi = static_cast<double>(j)*DY;
      for (int i=1; i<=NX-1; ++i) {
         const double xi = static_cast<double>(i)*DX;

         const int l = ij2l(i,j,NX);
 
         const double ref = sinh(M_PI*yi)*sin(M_PI*xi)/sinh(M_PI);
         err = max(err,abs(ref-sol[l]));
         
      }
   }
   return err;
}
