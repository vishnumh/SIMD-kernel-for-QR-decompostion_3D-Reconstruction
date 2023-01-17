#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>



#define FMADD(dest, src1, src2) \
  __asm__ __volatile__ (      \
  "vfmadd231pd %[rsrc1], %[rsrc2], %[rdest]\n"  \
    : [rdest] "+x"(dest)     \
    : [rsrc1] "x" (src1) , [rsrc2] "x"(src2));
    

//timing routine for reading the time stamp counter
static __inline__ unsigned long long rdtsc(void) {
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}


void naive
(
  int              num_iter,
  int              m,
  int              n,
  double* restrict A,      
  double* restrict F,
  double* restrict Q
){
 
    int i, j, k;
    double temp11, temp21, temp12, temp22, temp13, temp23, temp14, temp24, ulen1, ulen2, ulen3, ulen4, vnorm1, vnorm2, vnorm3, vnorm4;

    double * a;
    double * ac;
    double * Qc;
    posix_memalign((void**) &a, 64, num_iter * 4 * n * n * sizeof(double));
    posix_memalign((void**) &Qc, 64, num_iter *4 * n * n * sizeof(double));
    posix_memalign((void**) &ac, 64, num_iter *4 * n * n * sizeof(double));

      for(int iter = 0 ; iter < num_iter; ++iter ){
            for (int j = 0; j < n; j++){
                for (int p = 0; p < n; p++){
                  
                  for (int k = 0; k < num_iter *4; k++){
                  a[num_iter*4*n*j + num_iter*4*p + k + 4*iter] = 0;
                  for (int i = 0; i < m; i++){
                        
                a[num_iter*4*n*j + num_iter*4*p + k + 4*iter] += A[num_iter*4*n*i + num_iter*4*j + k + 4*iter] * A[num_iter*4*n*i + num_iter*4*p + k + 4*iter]; 
                
                }
              
              }
              
            }
        
          }
    
    
   m = n;
    
    // Loop across i starts
    for(i = 0; i < n; i ++){
      
      ulen1 = 0;
      ulen2 = 0;
      ulen3 = 0;
      ulen4 = 0;
      
      for(j = i; j < m; j ++){
        
        ulen1 += a[num_iter*4*n*j + num_iter*4*i + 0 + 4*iter]*a[num_iter*4*n*j + num_iter*4*i + 0 + 4*iter];
        ulen2 += a[num_iter*4*n*j + num_iter*4*i + 1 + 4*iter]*a[num_iter*4*n*j + num_iter*4*i + 1 + 4*iter];
        ulen3 += a[num_iter*4*n*j + num_iter*4*i + 2 + 4*iter]*a[num_iter*4*n*j + num_iter*4*i + 2 + 4*iter];
        ulen4 += a[num_iter*4*n*j + num_iter*4*i + 3 + 4*iter]*a[num_iter*4*n*j + num_iter*4*i + 3 + 4*iter];
      }
      ulen1 = sqrt(ulen1);
      ulen2 = sqrt(ulen2);
      ulen3 = sqrt(ulen3);
      ulen4 = sqrt(ulen4);
      
      // Initial vnorm
      if(a[num_iter*4*n*i + num_iter*4*i + 0 + 4*iter] < 0)
        ulen1 = ulen1;
      else
        ulen1 = -ulen1;
        
      if(a[num_iter*4*n*i + num_iter*4*i + 1 + 4*iter] < 0)
        ulen2 = ulen2;
      else
        ulen2 = -ulen2;
        
      if(a[num_iter*4*n*i + num_iter*4*i + 2 + 4*iter] < 0)
        ulen3 = ulen3;
      else
        ulen3 = -ulen3;
        
      if(a[num_iter*4*n*i + num_iter*4*i + 3 + 4*iter] < 0)
        ulen4 = ulen4;
      else
        ulen4 = -ulen4;

      
      // Modified vnorm
      vnorm1 = 0;
      vnorm2 = 0;
      vnorm3 = 0;
      vnorm4 = 0;
      
      for(j = i; j < m; j ++){

        if(j == i){
          vnorm1 += (a[num_iter*4*n*j + num_iter*4*i + 0 + 4*iter] + ulen1)*(a[num_iter*4*n*j + num_iter*4*i+ 0 + 4*iter] + ulen1);
          vnorm2 += (a[num_iter*4*n*j + num_iter*4*i + 1 +4*iter] + ulen2)*(a[num_iter*4*n*j + num_iter*4*i+ 1 +4*iter] + ulen2);
          vnorm3 += (a[num_iter*4*n*j + num_iter*4*i + 2 + 4*iter] + ulen3)*(a[num_iter*4*n*j + num_iter*4*i+ 2 + 4*iter] + ulen3);
          vnorm4 += (a[num_iter*4*n*j + num_iter*4*i + 3 +4*iter] + ulen4)*(a[num_iter*4*n*j + num_iter*4*i+ 3 +4*iter] + ulen4);
        
        }
        else{
          vnorm1 += (a[num_iter*4*n*j + num_iter*4*i + 0 + 4*iter])*(a[num_iter*4*n*j + num_iter*4*i + 0 + 4*iter]);
          vnorm2 += (a[num_iter*4*n*j + num_iter*4*i + 1 + 4*iter])*(a[num_iter*4*n*j + num_iter*4*i + 1 + 4*iter]);
          vnorm3 += (a[num_iter*4*n*j + num_iter*4*i + 2 + 4*iter])*(a[num_iter*4*n*j + num_iter*4*i + 2 + 4*iter]);
          vnorm4 += (a[num_iter*4*n*j + num_iter*4*i + 3 + 4*iter])*(a[num_iter*4*n*j + num_iter*4*i + 3 +4*iter]);
        }
      }
      
      vnorm1 = sqrt(vnorm1);
      vnorm2 = sqrt(vnorm2);
      vnorm3 = sqrt(vnorm3);
      vnorm4 = sqrt(vnorm4);
      
      if(vnorm1 < 1e-18){
        i = i + 1;
        if(i >= n){
          break;
        }
      }
      
      // Copy value for matrix multiplication
      for(j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
          
          ac[num_iter*4*j*m + num_iter*4*k + 0 +4*iter] = a[num_iter*4*j*m + num_iter*4*k + 0 +4*iter];
          ac[num_iter*4*j*m + num_iter*4*k + 1 +4*iter] = a[num_iter*4*j*m + num_iter*4*k + 1 +4*iter];
          ac[num_iter*4*j*m + num_iter*4*k + 2 +4*iter] = a[num_iter*4*j*m + num_iter*4*k + 2 +4*iter];
          ac[num_iter*4*j*m + num_iter*4*k + 3 +4*iter] = a[num_iter*4*j*m + num_iter*4*k + 3 +4*iter];
          
//          if(i == 0){
//            if(j == k){
//              Qc[num_iter*4*j*m + num_iter*4*k + 0+4*iter] = 1;
//              Qc[num_iter*4*j*m + num_iter*4*k + 1+4*iter] = 1;
//              Qc[num_iter*4*j*m + num_iter*4*k + 2+4*iter] = 1;
//              Qc[num_iter*4*j*m + num_iter*4*k + 3+4*iter] = 1;}
//            else{
//              Qc[num_iter*4*j*m + num_iter*4*k + 0+4*iter] = 0;
//              Qc[num_iter*4*j*m + num_iter*4*k + 1+4*iter] = 0;
//              Qc[num_iter*4*j*m + num_iter*4*k + 2+4*iter] = 0;
//              Qc[num_iter*4*j*m + num_iter*4*k + 3+4*iter] = 0;}
//          }
//          else{
//          Qc[num_iter*4*j*m + num_iter*4*k + 0+4*iter] = Q[num_iter*4*j*m + num_iter*4*k + 0+ 4*iter];
//          Qc[num_iter*4*j*m + num_iter*4*k + 1+4*iter] = Q[num_iter*4*j*m + num_iter*4*k + 1+ 4*iter];
//          Qc[num_iter*4*j*m + num_iter*4*k + 2+4*iter] = Q[num_iter*4*j*m + num_iter*4*k + 2+ 4*iter];
//          Qc[num_iter*4*j*m + num_iter*4*k + 3+4*iter] = Q[num_iter*4*j*m + num_iter*4*k + 3+ 4*iter];}
        
        }
      }
      
//      // Q = Q * vvT (ignoring zero multiplications)
//      for(int l = i; l < m; l ++){
//        for(int j = 0; j < n; j++){
//        
//          Q[num_iter*4*l*m + num_iter*4*j + 0 + 4*iter] = 0;
//          Q[num_iter*4*l*m + num_iter*4*j + 1 + 4*iter] = 0;
//          Q[num_iter*4*l*m + num_iter*4*j + 2 + 4*iter] = 0;
//          Q[num_iter*4*l*m + num_iter*4*j + 3+ 4*iter] = 0;
//        
//          for(int k = i; k < m ; k++){
//          
//            if(l == i){
//            temp11 = ac[num_iter*4*n*l + num_iter*4*i + 0 + 4*iter] + ulen1;
//            temp12 = ac[num_iter*4*n*l + num_iter*4*i + 1 + 4*iter] + ulen2;
//            temp13 = ac[num_iter*4*n*l + num_iter*4*i + 2 + 4*iter] + ulen3;
//            temp14 = ac[num_iter*4*n*l + num_iter*4*i + 3 + 4*iter] + ulen4;}
//            
//            else{
//            temp11 = ac[num_iter*4*n*l + num_iter*4*i + 0 + 4*iter];
//            temp12 = ac[num_iter*4*n*l + num_iter*4*i + 1 + 4*iter];
//            temp13 = ac[num_iter*4*n*l + num_iter*4*i + 2 + 4*iter];
//            temp14 = ac[num_iter*4*n*l + num_iter*4*i + 3 + 4*iter];}
//            
//            
//            if(k == i){
//            temp21 = ac[num_iter*4*n*k + num_iter*4*i + 0 + 4*iter] + ulen1; 
//            temp22 = ac[num_iter*4*n*k + num_iter*4*i + 1 + 4*iter] + ulen2;
//            temp23 = ac[num_iter*4*n*k + num_iter*4*i + 2 + 4*iter] + ulen3;
//            temp24 = ac[num_iter*4*n*k + num_iter*4*i + 3 + 4*iter] + ulen4;}
//            
//            else{
//            temp21 = ac[num_iter*4*n*k + num_iter*4*i + 0 + 4*iter];
//            temp22 = ac[num_iter*4*n*k + num_iter*4*i + 1 + 4*iter];
//            temp23 = ac[num_iter*4*n*k + num_iter*4*i + 2 + 4*iter];
//            temp24 = ac[num_iter*4*n*k + num_iter*4*i + 3 + 4*iter];}
//          
//            if(l == k){
//              temp11 = 1 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
//              temp12 = 1 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
//              temp13 = 1 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
//              temp14 = 1 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
//                      
//            else{
//              temp11 = 0 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
//              temp12 = 0 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
//              temp13 = 0 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
//              temp14 = 0 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
//            
//            Q[num_iter*4*l*m + num_iter*4*j + 0 + 4*iter] += temp11 * Qc[num_iter*4*k*m + num_iter*4*j + 0 + 4*iter];
//            Q[num_iter*4*l*m + num_iter*4*j + 1 + 4*iter] += temp12 * Qc[num_iter*4*k*m + num_iter*4*j + 1 + 4*iter];
//            Q[num_iter*4*l*m + num_iter*4*j + 2 + 4*iter] += temp13 * Qc[num_iter*4*k*m + num_iter*4*j + 2 + 4*iter];
//            Q[num_iter*4*l*m + num_iter*4*j + 3 + 4*iter] += temp14 * Qc[num_iter*4*k*m + num_iter*4*j + 3 + 4*iter];          
//          }
//          
//               
//        }
//
//      }
      
      // R = R * vvT (ignoring zero multiplications)
      for(int l = i; l < m; l ++){
        for(int j = i; j < n; j++){
        
          a[num_iter*4*l*m + num_iter*4*j + 0 + 4*iter] = 0;
          a[num_iter*4*l*m + num_iter*4*j + 1 + 4*iter] = 0;
          a[num_iter*4*l*m + num_iter*4*j + 2 + 4*iter] = 0;
          a[num_iter*4*l*m + num_iter*4*j + 3 + 4*iter] = 0;
          
        
          for(int k = i; k < m ; k++){
          
            if(l == i){
            temp11 = ac[num_iter*4*n*l + num_iter*4*i + 0 + 4*iter] + ulen1;
            temp12 = ac[num_iter*4*n*l + num_iter*4*i + 1+ 4*iter] + ulen2;
            temp13 = ac[num_iter*4*n*l + num_iter*4*i + 2+ 4*iter] + ulen3;
            temp14 = ac[num_iter*4*n*l + num_iter*4*i + 3+ 4*iter] + ulen4;}
            
            else{
            temp11 = ac[num_iter*4*n*l + num_iter*4*i + 0+ 4*iter];
            temp12 = ac[num_iter*4*n*l + num_iter*4*i + 1+ 4*iter];
            temp13 = ac[num_iter*4*n*l + num_iter*4*i + 2+ 4*iter];
            temp14 = ac[num_iter*4*n*l + num_iter*4*i + 3+ 4*iter];}
            
            if(k == i){
            temp21 = ac[num_iter*4*n*k + num_iter*4*i + 0+ 4*iter] + ulen1;
            temp22 = ac[num_iter*4*n*k + num_iter*4*i + 1+ 4*iter] + ulen2;
            temp23 = ac[num_iter*4*n*k + num_iter*4*i + 2+ 4*iter] + ulen3;
            temp24 = ac[num_iter*4*n*k + num_iter*4*i + 3+ 4*iter] + ulen4;}
            
            else{
            temp21 = ac[num_iter*4*n*k + num_iter*4*i + 0+ 4*iter];
            temp22 = ac[num_iter*4*n*k + num_iter*4*i + 1+ 4*iter];
            temp23 = ac[num_iter*4*n*k + num_iter*4*i + 2+ 4*iter];
            temp24 = ac[num_iter*4*n*k + num_iter*4*i + 3+ 4*iter];}
                      
            if(l == k){
              temp11 = 1 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
              temp12 = 1 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
              temp13 = 1 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
              temp14 = 1 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
                        
            else{
            
              temp11 = 0 - 2 * (temp11/vnorm1) * (temp21/vnorm1);
              temp12 = 0 - 2 * (temp12/vnorm2) * (temp22/vnorm2);
              temp13 = 0 - 2 * (temp13/vnorm3) * (temp23/vnorm3);
              temp14 = 0 - 2 * (temp14/vnorm4) * (temp24/vnorm4);}
        
            a[num_iter*4*l*m + num_iter*4*j + 0+ 4*iter] += temp11 * ac[num_iter*4*k*m + num_iter*4*j + 0+ 4*iter];
            a[num_iter*4*l*m + num_iter*4*j + 1+ 4*iter] += temp12 * ac[num_iter*4*k*m + num_iter*4*j + 1+ 4*iter];
            a[num_iter*4*l*m + num_iter*4*j + 2+ 4*iter] += temp13 * ac[num_iter*4*k*m + num_iter*4*j + 2+ 4*iter];
            a[num_iter*4*l*m + num_iter*4*j + 3+ 4*iter] += temp14 * ac[num_iter*4*k*m + num_iter*4*j + 3+ 4*iter];
                    
          }
         }
        }
      }
    }
  }



void matrixmulkernel(
  int              num_iter,
  int              iter,
  int              k,
  int              n,
  int              m,
  double* restrict A,
  double* restrict At
){
  __m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
  __m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;
  //double* aptr = A;
//  ymm0 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 0 * n + 0 ] );  ymm1 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 0 * n + 4 ] );
//  ymm2 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k +1 * n + 0 ] );  ymm3 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 1 * n + 4 ] );
//  ymm4 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 2 * n + 0 ] );  ymm5 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 2 * n + 4 ] );
//  ymm6 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 3 * n + 0 ] );  ymm7 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 3 * n + 4 ] );
//  ymm8 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 4* n + 0 ] );  ymm9 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 4 * n + 4 ] );
//  ymm10 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 5 * n + 0 ] );  ymm11 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 5 * n + 4 ] );
//  
  ymm0 = _mm256_setzero_pd();  ymm1 = _mm256_setzero_pd();
  ymm2 = _mm256_setzero_pd();  ymm3 = _mm256_setzero_pd();
  ymm4 = _mm256_setzero_pd();  ymm5 = _mm256_setzero_pd();
  ymm6 = _mm256_setzero_pd();  ymm7 = _mm256_setzero_pd();
  ymm8 = _mm256_setzero_pd();  ymm9 = _mm256_setzero_pd();
  ymm10 = _mm256_setzero_pd();  ymm11 = _mm256_setzero_pd();
  
  int i = 0;
  
  for (int p = 0; p < m; ++p)
    {
      ymm15 = _mm256_loadu_pd(&A[4*m*n*iter+ m*n*k + p * n + 0 ]);
      ymm14 = _mm256_loadu_pd(&A[4*m*n*iter+ m*n*k + p * n + 4 ]);
  
      //-------------------------------------------
      ymm13 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 0] );
      ymm12 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 1] );
      
      ymm0 = _mm256_fmadd_pd(ymm13, ymm15, ymm0);
      ymm1 = _mm256_fmadd_pd(ymm13, ymm14, ymm1);
    
      
      //-------------------------------------------
      ymm2 = _mm256_fmadd_pd(ymm12, ymm15, ymm2);
      ymm3 = _mm256_fmadd_pd(ymm12, ymm14, ymm3);
    
      
      //-------------------------------------------
      ymm13 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 2] );
      ymm12 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 3] );
      
      ymm4 = _mm256_fmadd_pd(ymm13, ymm15, ymm4);
      ymm5 = _mm256_fmadd_pd(ymm13, ymm14, ymm5);
       
      
      //-------------------------------------------
      ymm6 = _mm256_fmadd_pd(ymm12, ymm15, ymm6);
      ymm7 = _mm256_fmadd_pd(ymm12, ymm14, ymm7);
       
      
      //-------------------------------------------
      ymm13 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 4] );
      ymm12 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 5] );
      
      ymm8 = _mm256_fmadd_pd(ymm13, ymm15, ymm8);
      ymm9 = _mm256_fmadd_pd(ymm13, ymm14, ymm9);
    
      
      ymm10 = _mm256_fmadd_pd(ymm12, ymm15, ymm10);
      ymm11 = _mm256_fmadd_pd(ymm12, ymm14, ymm11);
     
      
    }


    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  0 + n*0 , ymm0);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  4 + n*0 , ymm1);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  0 + n*1 , ymm2);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  4 + n*1 , ymm3);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  0 + n*2 , ymm4);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  4 + n*2 , ymm5);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  0 + n*3 , ymm6);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  4 + n*3 , ymm7);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  0 + n*4 , ymm8);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  4 + n*4 , ymm9);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  0 + n*5 , ymm10);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k +  4 + n*5 , ymm11);

//  ymm0 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 6 * n + 0 ] );  ymm1 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 6 * n + 4 ] );
//  ymm2 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 7 * n + 0 ] );  ymm3 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 7 * n + 4 ] );
//  ymm4 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 8 * n + 0 ] );  ymm5 = _mm256_loadu_pd(&At[4*n*n*iter+ n*n*k + 8 * n + 4 ] );
  ymm0 = _mm256_setzero_pd();  ymm1 = _mm256_setzero_pd();
  ymm2 = _mm256_setzero_pd();  ymm3 = _mm256_setzero_pd();
  ymm4 = _mm256_setzero_pd();  ymm5 = _mm256_setzero_pd();

  for (int p = 0; p < m; ++p)
    {
      ymm15 = _mm256_loadu_pd(&A[4*m*n*iter+ m*n*k + p * n + 0 ]);
      ymm14 = _mm256_loadu_pd(&A[4*m*n*iter+ m*n*k + p * n + 4 ]);
      //-------------------------------------------
      ymm13 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 6 ] );
      ymm12 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 7 ] );
      ymm11 = _mm256_broadcast_sd(&A[4*m*n*iter+ m*n*k + p * n + 8 ] );
      ymm0 = _mm256_fmadd_pd(ymm13, ymm15, ymm0);
      ymm1 = _mm256_fmadd_pd(ymm13, ymm14, ymm1);
      ymm2 = _mm256_fmadd_pd(ymm12, ymm15, ymm2);
      ymm3 = _mm256_fmadd_pd(ymm12, ymm14, ymm3);
      ymm4 = _mm256_fmadd_pd(ymm11, ymm15, ymm4);
      ymm5 = _mm256_fmadd_pd(ymm11, ymm14, ymm5);
    }
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k + 0 + n*6 , ymm0);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k + 4 + n*6 , ymm1);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k + 0 + n*7 , ymm2);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k + 4 + n*7 , ymm3);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k + 0 + n*8 , ymm4);
    _mm256_storeu_pd(At + 4*n*n*iter+ n*n*k + 4 + n*8 , ymm5);
  }



void kernel
(
  int              num_iter,
  int              iter,
  int              m,
  int              n,      
  double* restrict a,
  double* restrict ac
  
){
    
   
    unsigned long long t1, t2, t3, t4 , tt, tb;
    //tt = 0; tb = 0;
    __m256d asimd1, asimd2, asimd3, asimd4, asimd5, asimd6, asimd7, asimd8, asimd9, asimd10, asimd11, asimd12, asimd13, asimd14, asimd15, asimd16;

    // 4*n*n*iter +  4*n*i + 4*j + k
    m = n;
    
    // Loop across i starts
    for(int i = 0; i < n; i ++){
    
      t1 = rdtsc();
      //printf("iteration i, %d\n", i);
      
      asimd11 = _mm256_setzero_pd();
      asimd1  = _mm256_load_pd(a + 4*n*n*iter +  4*n*0 + 4*i);
      asimd2  = _mm256_load_pd(a + 4*n*n*iter +  4*n*1 + 4*i);
      asimd3  = _mm256_load_pd(a + 4*n*n*iter +  4*n*2 + 4*i);
      asimd4  = _mm256_load_pd(a + 4*n*n*iter +  4*n*3 + 4*i);
      asimd5  = _mm256_load_pd(a + 4*n*n*iter +  4*n*4 + 4*i);
      
      if(i<1)
      FMADD(asimd11, asimd1, asimd1);
      if(i<2)
      FMADD(asimd11, asimd2, asimd2);
      if(i<3)
      FMADD(asimd11, asimd3, asimd3);
      if(i<4)
      FMADD(asimd11, asimd4, asimd4);
      if(i<5)
      FMADD(asimd11, asimd5, asimd5);
      
      
      asimd6  = _mm256_load_pd(a + 4*n*n*iter +  4*n*5 + 4*i);
      asimd7  = _mm256_load_pd(a + 4*n*n*iter +  4*n*6 + 4*i);
      asimd8  = _mm256_load_pd(a + 4*n*n*iter +  4*n*7 + 4*i);
      asimd9  = _mm256_load_pd(a + 4*n*n*iter +  4*n*8 + 4*i);
      
      if(i<6)
      FMADD(asimd11, asimd6, asimd6);
      if(i<7)
      FMADD(asimd11, asimd7, asimd7);
      if(i<8)
      FMADD(asimd11, asimd8, asimd8);
      if(i<9)
      FMADD(asimd11, asimd9, asimd9);


      
      asimd10 = _mm256_load_pd(a + 4*n*n*iter +  4*n*i + 4*i);
      asimd13 = _mm256_andnot_pd (asimd10, _mm256_set1_pd(-0.0));
      asimd12 = _mm256_sub_pd(asimd11, _mm256_mul_pd(asimd10, asimd10));
      asimd11 = _mm256_sqrt_pd(asimd11);
      asimd11 = _mm256_xor_pd(asimd11, asimd13);
      asimd10 = _mm256_add_pd(asimd10, asimd11);

      FMADD(asimd12, asimd10, asimd10);
      
      asimd12 = _mm256_add_pd(asimd12, _mm256_set1_pd(1.0*1e-24));
      asimd12 = _mm256_div_pd(_mm256_set1_pd(1.0), asimd12);
      
      t2 = rdtsc();
      // --------------------------------------------------------------------
      // Copy value for matrix multiplication
      for(int j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
          //num_iter*4*n*8 + num_iter*4*i + 4*iter
          
          ac[4*n*n*iter +  4*n*j + 4*k + 0] = a[4*n*n*iter +  4*n*j + 4*k + 0];
          ac[4*n*n*iter +  4*n*j + 4*k + 1] = a[4*n*n*iter +  4*n*j + 4*k + 1];
          ac[4*n*n*iter +  4*n*j + 4*k + 2] = a[4*n*n*iter +  4*n*j + 4*k + 2];
          ac[4*n*n*iter +  4*n*j + 4*k + 3] = a[4*n*n*iter +  4*n*j + 4*k + 3];
        
        }
      }
    //---------------------------------------------------------------------------
    t3 = rdtsc();
    for(int l = i; l < m; l ++){
            asimd16 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*l + 4*i);
              if(l == i){
                asimd16 = _mm256_add_pd(asimd16, asimd11);
              }
            for(int j = i; j < n; j++){
              asimd15 = _mm256_setzero_pd();
              if(i == 0){
              
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*0 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*1 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*2 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*i);    
              asimd5 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*i);    

              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11));
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              asimd10= _mm256_mul_pd(asimd16, asimd5);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              asimd10= _mm256_add_pd(asimd10, asimd10);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              asimd10= _mm256_mul_pd(asimd10, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 0)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd9);
              asimd10= _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd10);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*0 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*1 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*2 + 4*j));
              FMADD(asimd15, asimd9, _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*j));
              FMADD(asimd15, asimd10,_mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*j));
              
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              
              // 10 14  
              asimd6 = _mm256_mul_pd(asimd16, asimd1);
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd9);
              
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd9, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              
              }
              
              if(i == 1){
              
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*1 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*2 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11) );
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd9);
              
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*1 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*2 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*j));
              FMADD(asimd15, asimd9,_mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*j));
              
              
              
              // -----------------------------------------------------
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, asimd1);
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd9);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd9, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              
              }
              
              
              if(i == 2){
              
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*2 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11) );
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd8);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*2 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*j));
              FMADD(asimd15, asimd8,_mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*j));
              
              
              
              // -----------------------------------------------------
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, asimd1);
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd9);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd9, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              
              
              }
              
              
              if(i == 3){
              
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*i);   
              
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11));
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd7);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*3 + 4*j));
              FMADD(asimd15, asimd7,_mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*j));
              
              
              // -----------------------------------------------------
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, asimd1);
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd9);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd9, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              
              }
              
              
              if(i == 4){
              
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*i);    
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11));
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd6);
              FMADD(asimd15, asimd6,_mm256_load_pd(ac + 4*n*n*iter +  4*n*4 + 4*j));
              
              // -----------------------------------------------------
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, asimd1);
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd9);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd9, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              
              }
              
              if(i == 5){
               // -----------------------------------------------------
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);    
              asimd4 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11));
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              asimd9 = _mm256_mul_pd(asimd16, asimd4);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              asimd9 = _mm256_add_pd(asimd9, asimd9);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              asimd9 = _mm256_mul_pd(asimd9, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd8);
              asimd9 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd9);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*5 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd9, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              
              }
              
              if(i == 6){
              
               // -----------------------------------------------------
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);   
              asimd3 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11));
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              asimd8 = _mm256_mul_pd(asimd16, asimd3);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              asimd8 = _mm256_add_pd(asimd8, asimd8);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              asimd8 = _mm256_mul_pd(asimd8, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd7);
              asimd8 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd8);
              
             
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*6 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd8, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              }
              
              if(i == 7){
               // -----------------------------------------------------
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*i);    
              asimd2 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);   
              
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11));
              asimd7 = _mm256_mul_pd(asimd16, asimd2);
              
              
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd7 = _mm256_add_pd(asimd7, asimd7);
              
              
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd7 = _mm256_mul_pd(asimd7, asimd12);
              
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd6);
              asimd7 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd7);
              
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*7 + 4*j));
              FMADD(asimd15, asimd7, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              
              }
              
              if(i == 8){
              asimd1 = _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*i);    
              asimd6 = _mm256_mul_pd(asimd16, _mm256_add_pd(asimd1, asimd11));
              asimd6 = _mm256_add_pd(asimd6, asimd6);
              asimd6 = _mm256_mul_pd(asimd6, asimd12);
              asimd6 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd6);
              FMADD(asimd15, asimd6, _mm256_load_pd(ac + 4*n*n*iter +  4*n*8 + 4*j));
              }
          
          _mm256_store_pd(a + 4*n*n*iter + n*4*l + 4*j , asimd15); 
        }
      }
      t4 = rdtsc();
      tt += (t2-t1);
      tb += (t4-t3);
    
    } 
    
if(iter ==num_iter - 1)printf(" \n Top half Kernel Cyles = %lf \n",  ((double)(t2-t1)));

//if(iter ==num_iter - 1)printf(" \n Top half Kernel Cyles = %lf \n",  ((double)((tt)/(9*1.0*num_iter ))));
//if(iter ==num_iter - 1)printf(" \n Bottom half Kernel Cyles = %lf \n",  ((double)(tb/(9*1.0* num_iter))));
//if(iter ==num_iter - 1)printf("\n lf \n", num_iter);
    // Print out the results
//    if(iter == 0){
//      printf("R = \n");
//      for(int i = 0; i < m; i++) {
//          for(int j = 0; j < n; j++) {
//              printf("%9.6g ", a[4*n*n*iter + n*4*i + 4*j]);
//          }
//          printf("\n");
//      }
//      printf("\n");
//    }

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
  
){
    
   
   

    __m256d asimd1, asimd2, asimd3, asimd4, asimd5, asimd6, asimd7, asimd8, asimd9, asimd10, asimd11, asimd12, asimd13, asimd14, asimd15, asimd16;
       
    // Find A^T A
    for(int iter = 0 ; iter < num_iter; ++iter ){
      for (int j = 0; j < n; j++){
          for (int p = 0; p < n; p++){
            
            asimd1  = _mm256_load_pd(a + num_iter*4*n*j + num_iter*4*p + 4*iter);
            asimd2 = _mm256_load_pd(A + num_iter*4*n*0 + num_iter*4*j + 4*iter);
            asimd3 = _mm256_load_pd(A + num_iter*4*n*0 + num_iter*4*p + 4*iter);
            asimd4 = _mm256_load_pd(A + num_iter*4*n*1 + num_iter*4*j + 4*iter);
            asimd5 = _mm256_load_pd(A + num_iter*4*n*1 + num_iter*4*p + 4*iter);
            asimd6 = _mm256_load_pd(A + num_iter*4*n*2 + num_iter*4*j + 4*iter);
            asimd7 = _mm256_load_pd(A + num_iter*4*n*2 + num_iter*4*p + 4*iter);
            FMADD(asimd1, asimd2, asimd3);
            FMADD(asimd1, asimd4, asimd5);
            FMADD(asimd1, asimd6, asimd7);
            
            asimd8 = _mm256_load_pd(A + num_iter*4*n*3 + num_iter*4*j + 4*iter);
            asimd9 = _mm256_load_pd(A + num_iter*4*n*3 + num_iter*4*p + 4*iter);
            asimd10 = _mm256_load_pd(A + num_iter*4*n*4 + num_iter*4*j + 4*iter);
            asimd11 = _mm256_load_pd(A + num_iter*4*n*4 + num_iter*4*p + 4*iter);
            asimd12 = _mm256_load_pd(A + num_iter*4*n*5 + num_iter*4*j + 4*iter);
            asimd13 = _mm256_load_pd(A + num_iter*4*n*5 + num_iter*4*p + 4*iter);
            asimd14 = _mm256_load_pd(A + num_iter*4*n*6 + num_iter*4*j + 4*iter);
            asimd15 = _mm256_load_pd(A + num_iter*4*n*6 + num_iter*4*p + 4*iter);

            FMADD(asimd1, asimd8, asimd9);
            FMADD(asimd1, asimd10, asimd11);
            FMADD(asimd1, asimd12, asimd13);
            
            asimd2 = _mm256_load_pd(A + num_iter*4*n*7 + num_iter*4*j + 4*iter);
            asimd3 = _mm256_load_pd(A + num_iter*4*n*7 + num_iter*4*p + 4*iter);
            FMADD(asimd1, asimd14, asimd15);
            FMADD(asimd1, asimd2, asimd3);

            
            _mm256_store_pd(a + num_iter*4*n*j + num_iter*4*p + 4*iter, asimd1);
            

        }
      }
      

   m = n;
    
    // Loop across i starts
    for(int i = 0; i < n; i ++){
      
      
      asimd11 = _mm256_setzero_pd();
      
      asimd1  = _mm256_load_pd(a + num_iter*4*n*0 + num_iter*4*i + 4*iter);
      asimd2  = _mm256_load_pd(a + num_iter*4*n*1 + num_iter*4*i + 4*iter);
      asimd3  = _mm256_load_pd(a + num_iter*4*n*2 + num_iter*4*i + 4*iter);
      asimd4  = _mm256_load_pd(a + num_iter*4*n*3 + num_iter*4*i + 4*iter);
      
      if(i<1)
      FMADD(asimd11, asimd1, asimd1);
      if(i<2)
      FMADD(asimd11, asimd2, asimd2);
      if(i<3)
      FMADD(asimd11, asimd3, asimd3);
      if(i<4)
      FMADD(asimd11, asimd4, asimd4);
      
      asimd5  = _mm256_load_pd(a + num_iter*4*n*4 + num_iter*4*i + 4*iter);
      asimd6  = _mm256_load_pd(a + num_iter*4*n*5 + num_iter*4*i + 4*iter);
      asimd7  = _mm256_load_pd(a + num_iter*4*n*6 + num_iter*4*i + 4*iter);
      asimd8  = _mm256_load_pd(a + num_iter*4*n*7 + num_iter*4*i + 4*iter);
      asimd9  = _mm256_load_pd(a + num_iter*4*n*8 + num_iter*4*i + 4*iter);
      
      if(i<5)
      FMADD(asimd11, asimd5, asimd5);
      if(i<6)
      FMADD(asimd11, asimd6, asimd6);
      if(i<7)
      FMADD(asimd11, asimd7, asimd7);
      if(i<8)
      FMADD(asimd11, asimd8, asimd8);
      if(i<9)
      FMADD(asimd11, asimd9, asimd9);


      asimd11 = _mm256_sqrt_pd(asimd11);
      asimd10 = _mm256_load_pd(a + num_iter*4*n*i + num_iter*4*i + 4*iter);
      //bit operations - and not instruction replace multiplication
      asimd12 = _mm256_mul_pd ( _mm256_and_pd(asimd10, _mm256_set1_pd(-0.0)),  _mm256_set1_pd(-1.0));
      asimd11 = _mm256_xor_pd(asimd11, asimd12);
      


      asimd12 = _mm256_setzero_pd();
      
      if(i == 0){
      asimd1 = _mm256_add_pd(asimd1, asimd11);
      FMADD(asimd12, asimd1, asimd1);
      FMADD(asimd12, asimd2, asimd2);
      FMADD(asimd12, asimd3, asimd3);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 1){
      asimd2 = _mm256_add_pd(asimd2, asimd11);
      FMADD(asimd12, asimd2, asimd2);
      FMADD(asimd12, asimd3, asimd3);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 2){
      asimd3 = _mm256_add_pd(asimd3, asimd11);
      FMADD(asimd12, asimd3, asimd3);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 3){
      asimd4 = _mm256_add_pd(asimd4, asimd11);
      FMADD(asimd12, asimd4, asimd4);     
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 4){
      asimd5 = _mm256_add_pd(asimd5, asimd11);    
      FMADD(asimd12, asimd5, asimd5);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 5){
      asimd6 = _mm256_add_pd(asimd6, asimd11);
      FMADD(asimd12, asimd6, asimd6);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 6){
      asimd7 = _mm256_add_pd(asimd7, asimd11);
      FMADD(asimd12, asimd7, asimd7);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 7){
      asimd8 = _mm256_add_pd(asimd8, asimd11);
      FMADD(asimd12, asimd8, asimd8);
      FMADD(asimd12, asimd9, asimd9);}
      else if(i == 8){
      asimd9 = _mm256_add_pd(asimd9, asimd11);
      FMADD(asimd12, asimd9, asimd9);}

      
      asimd12 = _mm256_add_pd(asimd12, _mm256_set1_pd(1.0*1e-24));
      
      
      // Copy value for matrix multiplication
      for(int j = 0; j < m; j ++){
        for(int k = 0; k < m; k++){
          
          ac[num_iter*4*j*m + num_iter*4*k + 0 + 4*iter] = a[num_iter*4*j*m + num_iter*4*k + 0 + 4*iter];
          ac[num_iter*4*j*m + num_iter*4*k + 1 + 4*iter] = a[num_iter*4*j*m + num_iter*4*k + 1 + 4*iter];
          ac[num_iter*4*j*m + num_iter*4*k + 2 + 4*iter] = a[num_iter*4*j*m + num_iter*4*k + 2 + 4*iter];
          ac[num_iter*4*j*m + num_iter*4*k + 3 + 4*iter] = a[num_iter*4*j*m + num_iter*4*k + 3 + 4*iter];
        
        }
      }
      

    for(int l = i; l < m; l ++){
            asimd16 = _mm256_load_pd(ac + num_iter*4*n*l + num_iter*4*i + 4*iter);
              if(l == i){
                asimd16 = _mm256_add_pd(asimd16, asimd11);
              }
            for(int j = i; j < n; j++){

              
              
              
              
            
              if(i == 0){
              
              asimd15 = _mm256_load_pd(ac + num_iter*4*0*m + num_iter*4*j + 4*iter);    
              asimd13 = _mm256_mul_pd(asimd16, asimd1);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 0)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*1*m + num_iter*4*j + 4*iter); 
              asimd13 = _mm256_mul_pd(asimd16, asimd2);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*2*m + num_iter*4*j + 4*iter);   
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*3*m + num_iter*4*j + 4*iter);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*4*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*5*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*6*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*7*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              }
              
              else if(i == 1){
              asimd15 = _mm256_load_pd(ac + num_iter*4*1*m + num_iter*4*j + 4*iter);  
                   
              asimd13 = _mm256_mul_pd(asimd16, asimd2);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 1)), asimd13);
              
              asimd15 = _mm256_mul_pd(asimd13, asimd15);

              asimd14 = _mm256_load_pd(ac + num_iter*4*2*m + num_iter*4*j + 4*iter);   
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*3*m + num_iter*4*j + 4*iter);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*4*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*5*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*6*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*7*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 2){
              asimd15 = _mm256_load_pd(ac + num_iter*4*2*m + num_iter*4*j + 4*iter);    
              asimd13 = _mm256_mul_pd(asimd16, asimd3);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 2)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*3*m + num_iter*4*j + 4*iter);   
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*4*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*5*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*6*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*7*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 3){
              asimd15 = _mm256_load_pd(ac + num_iter*4*3*m + num_iter*4*j + 4*iter);    
              asimd13 = _mm256_mul_pd(asimd16, asimd4);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 3)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*4*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*5*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*6*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*7*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              
              else if(i == 4){
              asimd15 = _mm256_load_pd(ac + num_iter*4*4*m + num_iter*4*j + 4*iter);    
              asimd13 = _mm256_mul_pd(asimd16, asimd5);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 4)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*5*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*6*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*7*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 5){
              asimd15 = _mm256_load_pd(ac + num_iter*4*5*m + num_iter*4*j + 4*iter);                
              asimd13 = _mm256_mul_pd(asimd16, asimd6);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 5)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*6*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*7*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 6){
              asimd15 = _mm256_load_pd(ac + num_iter*4*6*m + num_iter*4*j + 4*iter);    
              asimd13 = _mm256_mul_pd(asimd16, asimd7);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 6)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*7*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              FMADD(asimd15, asimd13, asimd14);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 7){
              asimd15 = _mm256_load_pd(ac + num_iter*4*7*m + 4*j + num_iter*4*iter);    
              asimd13 = _mm256_mul_pd(asimd16, asimd8);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 7)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);
              
              asimd14 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);  
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              FMADD(asimd15, asimd13, asimd14);}
              
              else if(i == 8){
              asimd15 = _mm256_load_pd(ac + num_iter*4*8*m + num_iter*4*j + 4*iter);    
              asimd13 = _mm256_mul_pd(asimd16, asimd9);
              asimd13 = _mm256_mul_pd(asimd13, _mm256_set1_pd(2.0));
              asimd13 = _mm256_div_pd(asimd13, asimd12);
              asimd13 = _mm256_sub_pd(_mm256_set1_pd(1.0*(l == 8)), asimd13);
              asimd15 = _mm256_mul_pd(asimd13, asimd15);}
              
          
          _mm256_store_pd(a + num_iter*4*l*m + num_iter*4*j + 4*iter, asimd15);    
        }

      }}}
}
