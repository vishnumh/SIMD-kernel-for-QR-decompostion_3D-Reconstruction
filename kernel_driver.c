#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "immintrin.h"

#define MAX_FREQ 3.2
#define BASE_FREQ 2.4

//timing routine for reading the time stamp counter
static __inline__ unsigned long long rdtsc(void) {
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

void old_kernel
(
  int              num_iter,
  int              m,
  int              n,
  double* restrict A,
  double* restrict F,
  double* restrict Q,
  double* restrict Qc,
  double* restrict a,
  double* restrict ac
 );
 
void matrixmulkernel(
  int              num_iter,
  int              iter,
  int              k,
  int              n,
  int              m,
  double* restrict A,      
  double* restrict a
  
);

void kernel
(
  int              num_iter,
  int              iter,
  int              m,
  int              n,      
  double* restrict a,
  double* restrict ac
  
);

void naive
(
  int              num_iter,
  int              m,
  int              n,
  double* restrict A,
  double* restrict F,
  double* restrict Q

 );
 
 
int main(int argc, char** argv){
  //TODO: Change this according to your calculations for the size of the kernel
  int m = 8;  //m is the number of rows of A
  int n = 9;  //n is the number of columns of A
  int num_iter = atoi(argv[1]);;
  double *Q;
  double *F;
  double *A;
  double *At;
  double *Qc;
  double *a;
  double *ac;


  unsigned long long t0, t1, t2, t3, t4, t5, t6, t7, tm, ts, te;
  tm = 0;
  
  printf("num_iter = %d \n", num_iter);
  
// Dataset loading code
//  int ind = 0, in = 0; 
//  double myvariable=0.0;
//  char ch;
//  FILE* ptr1;
//  ptr1 = fopen("./toybus_corresp_raw.txt", "r");
//  double pts1[239][2];
//  double pts2[239][2];
//  if (NULL == ptr1) {
//        printf("file can't be opened \n");
//    }
//  do {
//        if(ind%4 == 0)
//          //printf("\n");
//          in = 0;
//        fscanf(ptr1,"%lf",&myvariable);
//        //printf("%lf ",myvariable);
//        if((in == 0) ||(in == 1))
//          pts1[(int)ind/4][in] = myvariable;
//        if((in == 2) ||(in == 3))
//          pts2[(int)ind/4][in-2] = myvariable;
//          
//        // Checking if character is not EOF.
//        // If it is EOF stop eading.
//        ind++ ;
//        in ++;
//    } while (ind < 956);
//// 
////    // Closing the file
//    fclose(ptr1);  

  
  
  
  
  // Sample data
  double Am[4][8][9] = 
    {
    {
        {4.51249166e+00, 1.57616189e-01, -2.21716994e+00, -1.09471621e+00, -3.82371892e-02, 5.37878417e-01, -2.03524844e+00, -7.10889077e-02, 1.00000000e+00},
        {4.55803664e+00, -2.34418662e-01, 1.95273973e+00, -1.77873502e+00, 9.14798888e-02, -7.62040068e-01, 2.33417520e+00, -1.20046035e-01, 1.00000000e+00},
        {1.77002133e-01, 6.52226503e-03, -3.41981902e-01, -7.60576940e-02, -2.80261277e-03,  1.46949387e-01, -5.17577485e-01, -1.90719596e-02, 1.00000000e+00},
        {4.64699012e-01, -2.03392431e-01, 5.67007552e-01, -5.09929101e-01, 2.23189025e-01, -6.22195536e-01, 8.19564062e-01, -3.58712032e-01, 1.00000000e+00},
        {5.04420718e-01, -3.66671837e-01,  5.32046419e-01, -9.48466587e-01, 6.89456188e-01, -1.00041143e+00,  9.48076522e-01, -6.89172643e-01, 1.00000000e+00},
        {3.19842188e-02,  4.01381137e-02, -1.76711092e-01, 3.32296863e-02, 4.17010943e-02, -1.83592233e-01, -1.80997233e-01, -2.27139752e-01, 1.00000000e+00},
        {4.07584825e+00, -1.33173690e+00, -1.90251974e+00, -3.07245014e+00, 1.00388802e+00,  1.43415473e+00, -2.14234216e+00,  6.99985852e-01, 1.00000000e+00},
        {4.53285432e-02,  1.66627095e-01, -2.05315656e-01,  1.61924150e-01, 5.95230924e-01, -7.33435504e-01, -2.20774899e-01, -8.11565462e-01, 1.00000000e+00}
   },
   {
        {4.51249166e+00, 1.57616189e-01, -2.21716994e+00, -1.09471621e+00, -3.82371892e-02, 5.37878417e-01, -2.03524844e+00, -7.10889077e-02, 1.00000000e+00},
        {4.55803664e+00, -2.34418662e-01, 1.95273973e+00, -1.77873502e+00, 9.14798888e-02, -7.62040068e-01, 2.33417520e+00, -1.20046035e-01, 1.00000000e+00},
        {1.77002133e-01, 6.52226503e-03, -3.41981902e-01, -7.60576940e-02, -2.80261277e-03,  1.46949387e-01, -5.17577485e-01, -1.90719596e-02, 1.00000000e+00},
        {4.64699012e-01, -2.03392431e-01, 5.67007552e-01, -5.09929101e-01, 2.23189025e-01, -6.22195536e-01, 8.19564062e-01, -3.58712032e-01, 1.00000000e+00},
        {5.04420718e-01, -3.66671837e-01,  5.32046419e-01, -9.48466587e-01, 6.89456188e-01, -1.00041143e+00,  9.48076522e-01, -6.89172643e-01, 1.00000000e+00},
        {3.19842188e-02,  4.01381137e-02, -1.76711092e-01, 3.32296863e-02, 4.17010943e-02, -1.83592233e-01, -1.80997233e-01, -2.27139752e-01, 1.00000000e+00},
        {4.07584825e+00, -1.33173690e+00, -1.90251974e+00, -3.07245014e+00, 1.00388802e+00,  1.43415473e+00, -2.14234216e+00,  6.99985852e-01, 1.00000000e+00},
        {4.53285432e-02,  1.66627095e-01, -2.05315656e-01,  1.61924150e-01, 5.95230924e-01, -7.33435504e-01, -2.20774899e-01, -8.11565462e-01, 1.00000000e+00}
   },
   {
        {4.51249166e+00, 1.57616189e-01, -2.21716994e+00, -1.09471621e+00, -3.82371892e-02, 5.37878417e-01, -2.03524844e+00, -7.10889077e-02, 1.00000000e+00},
        {4.55803664e+00, -2.34418662e-01, 1.95273973e+00, -1.77873502e+00, 9.14798888e-02, -7.62040068e-01, 2.33417520e+00, -1.20046035e-01, 1.00000000e+00},
        {1.77002133e-01, 6.52226503e-03, -3.41981902e-01, -7.60576940e-02, -2.80261277e-03,  1.46949387e-01, -5.17577485e-01, -1.90719596e-02, 1.00000000e+00},
        {4.64699012e-01, -2.03392431e-01, 5.67007552e-01, -5.09929101e-01, 2.23189025e-01, -6.22195536e-01, 8.19564062e-01, -3.58712032e-01, 1.00000000e+00},
        {5.04420718e-01, -3.66671837e-01,  5.32046419e-01, -9.48466587e-01, 6.89456188e-01, -1.00041143e+00,  9.48076522e-01, -6.89172643e-01, 1.00000000e+00},
        {3.19842188e-02,  4.01381137e-02, -1.76711092e-01, 3.32296863e-02, 4.17010943e-02, -1.83592233e-01, -1.80997233e-01, -2.27139752e-01, 1.00000000e+00},
        {4.07584825e+00, -1.33173690e+00, -1.90251974e+00, -3.07245014e+00, 1.00388802e+00,  1.43415473e+00, -2.14234216e+00,  6.99985852e-01, 1.00000000e+00},
        {4.53285432e-02,  1.66627095e-01, -2.05315656e-01,  1.61924150e-01, 5.95230924e-01, -7.33435504e-01, -2.20774899e-01, -8.11565462e-01, 1.00000000e+00}
   },
   {
        {4.51249166e+00, 1.57616189e-01, -2.21716994e+00, -1.09471621e+00, -3.82371892e-02, 5.37878417e-01, -2.03524844e+00, -7.10889077e-02, 1.00000000e+00},
        {4.55803664e+00, -2.34418662e-01, 1.95273973e+00, -1.77873502e+00, 9.14798888e-02, -7.62040068e-01, 2.33417520e+00, -1.20046035e-01, 1.00000000e+00},
        {1.77002133e-01, 6.52226503e-03, -3.41981902e-01, -7.60576940e-02, -2.80261277e-03,  1.46949387e-01, -5.17577485e-01, -1.90719596e-02, 1.00000000e+00},
        {4.64699012e-01, -2.03392431e-01, 5.67007552e-01, -5.09929101e-01, 2.23189025e-01, -6.22195536e-01, 8.19564062e-01, -3.58712032e-01, 1.00000000e+00},
        {5.04420718e-01, -3.66671837e-01,  5.32046419e-01, -9.48466587e-01, 6.89456188e-01, -1.00041143e+00,  9.48076522e-01, -6.89172643e-01, 1.00000000e+00},
        {3.19842188e-02,  4.01381137e-02, -1.76711092e-01, 3.32296863e-02, 4.17010943e-02, -1.83592233e-01, -1.80997233e-01, -2.27139752e-01, 1.00000000e+00},
        {4.07584825e+00, -1.33173690e+00, -1.90251974e+00, -3.07245014e+00, 1.00388802e+00,  1.43415473e+00, -2.14234216e+00,  6.99985852e-01, 1.00000000e+00},
        {4.53285432e-02,  1.66627095e-01, -2.05315656e-01,  1.61924150e-01, 5.95230924e-01, -7.33435504e-01, -2.20774899e-01, -8.11565462e-01, 1.00000000e+00}
   }
     };
     
    //create memory aligned buffers

  posix_memalign((void**) &Q, 64, num_iter * 4* n * n * sizeof(double));
  posix_memalign((void**) &F, 64, num_iter * 4* n * sizeof(double));
  posix_memalign((void**) &A, 64, num_iter * 4 * m * n * sizeof(double));
  posix_memalign((void**) &At, 64, num_iter * 4 * n * n * sizeof(double));
  
  posix_memalign((void**) &Qc, 64, num_iter * 4* n * n * sizeof(double));
  posix_memalign((void**) &a, 64, num_iter * 4* n * n * sizeof(double));
  posix_memalign((void**) &ac, 64, num_iter * 4* n * n * sizeof(double));
   
   for (int iter = 0; iter < num_iter; iter++){   
    for (int k = 0; k <  4; k++){
      for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
        
        A[4*m*n*iter+ m*n*k +  n*i + j ] = Am[k][i][j] + k*1e-105;
          }
        }
      }
    }
    
    
  for (int iter = 0; iter < num_iter; iter++){  
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        for (int k = 0; k < 4; k++){
        
        if(i == j){
        Q[4*n*n*iter +  4*n*i + 4*j + k] = 1.0;}
        
        else{
        Q[4*n*n*iter +  4*n*i + 4*j + k] = 0.0;}
        
        a[4*n*n*iter +  4*n*i + 4*j + k] = 0.0;
        ac[4*n*n*iter +  4*n*i + 4*j + k] = 0.0;
        
        }
      }
    }
   }   


  // initialize At with 0.0
  for (int iter = 0; iter < num_iter; iter++){  
    for (int k = 0; k <  4; k++){
        for (int i = 0; i < n; i++){
          for (int j = 0; j < n; j++){
        
            At[ 4*n*n*iter+ n*n*k + i * n + j ] = 0.0;
          }
        }
    }
   }
    
   
   

// Matrix mul start ----------------------------------------------
tm = 0;
t0 = rdtsc();
//#pragma omp parallel for
for (int iter = 0; iter < num_iter; iter++){  
  for (int k = 0; k <  4; k++){
      ts = rdtsc();
      matrixmulkernel(num_iter, iter, k,  n, m, A, At);
      te = rdtsc();
      tm += te - ts;
      
//    scalar part
      for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
          At[ 4*n*n*iter+ n*n*k + n*i + 8 ] += A[4*m*n*iter+ m*n*k + i + j * n] * A[4*m*n*iter+ m*n*k + 8 + j * n] ;
          
        }
       }
    }
  }
t1 = rdtsc();
//printf(" \n Peek Mat mul Kernel Cyles = %lf \n",  (2*n*8*m)/(((double)tm/ (1.0 *num_iter*4))*(MAX_FREQ/BASE_FREQ)));
// Matrxi mul end -------------------------------------------------------------------------------  

// Changing Layout
for (int iter = 0; iter < num_iter; iter++){  
  for (int k = 0; k <  4; k++){     
    for (int i = 0; i < n; i++){
          for (int j = 0; j < n; j++){
        
          a[4*n*n*iter +  4*n*i + 4*j + k] = At[4*n*n*iter+ n*n*k + n*i + j ] ;
          //if(iter ==3 && k == 3)printf("%lf ", At[4*n*n*iter+ n*n*k + n*i + j ]);
        }
        //if(iter == 3 && k == 3) printf("\n");
       }
    }
  }
    free(At);   
    
  t2 = rdtsc();
//  #pragma omp parallel for  
  for (int iter = 0; iter < num_iter; iter++){   
    
    kernel(num_iter, iter, m, n, a, ac);
    
  }  
  t3 = rdtsc();
  
  t4 = rdtsc();
  //old_kernel(num_iter, m, n, A, F, Q, Qc, a, ac);    
  t5 = rdtsc();
  t6 = rdtsc();
  //naive(num_iter,m,n,A,F,Q);
  t7 = rdtsc();
  
  printf(" \n updated Mat Mul Kernel Cycle = %lf \n", ((double)tm/ (1.0 *num_iter*4)));
  printf(" \n updated QR Kernel Cycle = %lf \n", (double)(t3-t2));
  //printf(" \n old matmul and QR Kernel Cycle = %lf \n", (double)(t5-t4));
  //printf(" \n Naive Cycle = %lf \n", (double)(t7-t6));
  printf(" \n updated total Kernel Cycle = %lf \n", (double)(t1-t0) + (double)(t3-t2));
   
  return 0;
}