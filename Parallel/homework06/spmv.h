#ifndef SPMV_H
#define SPMV_H

#include "common.h"

void allocate_ell_gpu(unsigned int* col_ind, double* vals, int m, int n, 
                      int nnz, double* x, unsigned int** dev_col_ind, 
                      double** dev_vals);
void allocate_csr_gpu(unsigned int* row_ptr, unsigned int* col_ind, 
                      double* vals, int m, int n, int nnz, double* x, 
                      unsigned int** dev_row_ptr, unsigned int** dev_col_ind,
                      double** dev_vals);
void allocate_data_gpu(double* x, double* b, double** dx, double** db, 
                       double** drk, double** dpk, double** dap, double** z1,
                       double** z2, int m, int n);
void get_result_gpu(double* dev_b, double* b, int m);
template <class T>
void CopyData(
  T* input,
  unsigned int N,
  unsigned int dsize,
  T** d_in);
int cg_gpu_csr(unsigned int* row_ptr, unsigned int* col_ind, double* vals, 
               double *x, double* b, double* drk, double* dpk, double* dap, 
               double* z1, double* z2, int m, int n, int nnz, int max_iter, 
               double tol);
void free_gpu(unsigned int* drp, unsigned int* dci, unsigned int* dec, 
              double* dev, double* dx, double* db, double* drk, double* dpk, 
              double* dap, double* z1, double* z2);

void test_dot(int n, double* b, double *b1, double** d_b, double** d_b1,
	      double** z1, double** z2);

#endif

