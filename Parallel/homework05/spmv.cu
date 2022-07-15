#include <iostream>
#include <stdio.h>
#include <assert.h>

#include <helper_cuda.h>
#include <cooperative_groups.h>

#include "spmv.h"
#define THREADS_PER_BLOCK 256

template <class T>
__global__ void
spmv_kernel_ell(unsigned int* col_ind, T* vals, int m, int n, int nnz, 
                double* x, double* b)
{
  
  //shared memory for the reduction
  __shared__ double reduce[THREADS_PER_BLOCK];
 
  int row = blockIdx.x;
  int tid = threadIdx.x;
  
  
  if (tid == 0 && row == 0)
    printf("n %d\n\n\n\n\n\n\n", n);
  
  if (row < m)
    {

      int begin = row * n;
      int chunk = n/THREADS_PER_BLOCK;
      
     
      // have each tread sum chunk number of entries in matrix and store at there
      // index. There will be a remainder if THREADS_PER_BLOCK does not divide
      // n evenlly...just have the final thread do the remainder in serial.
      if (chunk == 0)
	{

	  if (tid < n)
	    reduce[tid] = vals[begin + tid] * x[col_ind[begin + tid]];
	  
	}
      
      // bug somewhere here
      else if (chunk == 1)
	{
	  if (tid < n)
	    reduce[tid] = vals[begin + tid] * x[col_ind[begin + tid]];

	  if (tid == THREADS_PER_BLOCK - 1)
	    {
	      for (int i = begin + tid + 1 ; i < begin + n ; i++)
		{
		  reduce[tid] = vals[i] * x[col_ind[i]];
		}
	    }
	}
      else
	{
	  reduce[tid] = 0;
	  int i;
	  for (i = begin + tid*chunk ; i < begin + (tid+1)*chunk ; i++)
	    {
	      reduce[tid] += vals[i] * x[col_ind[i]];
	    }
	  if (tid == (THREADS_PER_BLOCK - 1))
	    {
	      //printf("i is on row %d: %d %d\n", row ,i, n*(row+1));
	      while (i < n*(row+1))
		{
		  reduce[tid] += vals[i] * x[col_ind[i]];
		  i++;
		}
	    }
	  
	}
      
         
      // do sequential reduction (assuming number of threads is a power of 2
      for (int s = THREADS_PER_BLOCK/2 ; s > 0 ; s>>=1)
	{	     
	  if (tid < s)
	    reduce[tid] += reduce[tid + s];
	  __syncthreads();
	}
      
      if (tid == 0)
	b[row] = reduce[tid];
    }
  
  /*
  if (tid == 0 && row == 0)
    {
      printf("\n");
      for (int i = 0 ; i < m ; i++)
	printf("\t%f", b[i]);
      printf("\n");
    }
  */
  

}




void spmv_gpu_ell(unsigned int* col_ind, double* vals, int m, int n, int nnz, 
                  double* x, double* b)
{
    // timers
    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float elapsedTime;

    // GPU execution parameters
    unsigned int blocks = m; 
    unsigned int threads = THREADS_PER_BLOCK; 
    unsigned int shared = threads * sizeof(double);

    dim3 dimGrid(blocks, 1, 1);
    dim3 dimBlock(threads, 1, 1);

    checkCudaErrors(cudaEventRecord(start, 0));
    for(unsigned int i = 0; i < MAX_ITER; i++) {
        cudaDeviceSynchronize();
        spmv_kernel_ell<double><<<dimGrid, dimBlock, shared>>>(col_ind, vals, 
                                                               m, n, nnz, x, b);

    }
    checkCudaErrors(cudaEventRecord(stop, 0));
    checkCudaErrors(cudaEventSynchronize(stop));
    checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, stop));
    printf("  Exec time (per itr): %0.8f s\n", (elapsedTime / 1e3 / MAX_ITER));

}




void allocate_ell_gpu(unsigned int* col_ind, double* vals, int m, int n, int n_new, 
                      int nnz, double* x, unsigned int** dev_col_ind, 
                      double** dev_vals, double** dev_x, double** dev_b)
{
  
  CopyData(col_ind, m*n_new, sizeof(int), dev_col_ind);
  CopyData(vals, m*n_new, sizeof(double), dev_vals);
  CopyData(x, n, sizeof(double), dev_x);
  double* b = (double*) malloc(sizeof(double) * m);
  CopyData(b, m, sizeof(double), dev_b);
  free(b);
  
  

}

void allocate_csr_gpu(unsigned int* row_ptr, unsigned int* col_ind, 
                      double* vals, int m, int n, int nnz, double* x, 
                      unsigned int** dev_row_ptr, unsigned int** dev_col_ind,
                      double** dev_vals, double** dev_x, double** dev_b)
{

  CopyData(row_ptr, m+1, sizeof(int), dev_row_ptr);
  CopyData(col_ind, nnz, sizeof(int), dev_col_ind);
  CopyData(vals, nnz, sizeof(double), dev_vals);
  CopyData(x, n, sizeof(double), dev_x);
  double* b = (double*) malloc(sizeof(double) * m);
  CopyData(b, m, sizeof(double), dev_b);
  free(b);
  
    
}

void get_result_gpu(double* dev_b, double* b, int m)
{
    // timers
    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float elapsedTime;


    checkCudaErrors(cudaEventRecord(start, 0));
    checkCudaErrors(cudaMemcpy(b, dev_b, sizeof(double) * m, 
                               cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaEventRecord(stop, 0));
    checkCudaErrors(cudaEventSynchronize(stop));
    checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, stop));
    printf("  Pinned Host to Device bandwidth (GB/s): %f\n",
         (m * sizeof(double)) * 1e-6 / elapsedTime);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

template <class T>
void CopyData(
  T* input,
  unsigned int N,
  unsigned int dsize,
  T** d_in)
{
  // timers
  cudaEvent_t start;
  cudaEvent_t stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float elapsedTime;

  // Allocate pinned memory on host (for faster HtoD copy)
  T* h_in_pinned = NULL;
  checkCudaErrors(cudaMallocHost((void**) &h_in_pinned, N * dsize));
  assert(h_in_pinned);
  memcpy(h_in_pinned, input, N * dsize);

  // copy data
  checkCudaErrors(cudaMalloc((void**) d_in, N * dsize));
  checkCudaErrors(cudaEventRecord(start, 0));
  checkCudaErrors(cudaMemcpy(*d_in, h_in_pinned,
                             N * dsize, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaEventRecord(stop, 0));
  checkCudaErrors(cudaEventSynchronize(stop));
  checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, stop));
  printf("  Pinned Device to Host bandwidth (GB/s): %f\n",
         (N * dsize) * 1e-6 / elapsedTime);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}


template <class T>
__global__ void
spmv_kernel(unsigned int* row_ptr, unsigned int* col_ind, T* vals, 
            int m, int n, int nnz, double* x, double* b)
{

 
  // each block will do a reduction between a single row and the column vector x.
  int row = blockIdx.x;
  int tid = threadIdx.x;
  
  /*
  if (row == 0 && tid == 0)
    {
      printf("\n");
      for (int i = 0 ; i < m+1 ; i++)
	printf("\t%d", row_ptr[i]);
    }
  */

  //shared memory for the reduction
  __shared__ double reduce[THREADS_PER_BLOCK];

  // guard against block reading past rows of matrix
  // (unecessary if kernal is launched with m blocks)
  if (row < m)
    {
      
      // set reduction array to zero
      reduce[tid] = 0;
      
   

      int begin = row_ptr[row];
      int end = row_ptr[row+1];
      int length = end - begin;

     
      /*
      if (tid == 0)
	printf("row %d length: %d\n", row, length);
      */
      //fill reduction array with products
      // since reduction array is only length THREADS_PER_BLOCK
      // When we first put products in the array they must be summed down to fit inside.
      // This can be done by taking chunks of the product vector and having one
      // thread sum them all
      // these are stored in the thread id index
      // Once this is done a sequential reduction can be done.
      reduce[tid] = 0;
      int chunk = length/THREADS_PER_BLOCK;
      //if(row == 5 && tid == 0)
      //	printf("chunk is %d %d\n", chunk, length);
      // if we can fit everything in row directly into shared.
      if (chunk == 0)
	{
	  //printf("in 0 chunk\n");
	  if (tid < length)
	    reduce[tid] = vals[begin + tid] * x[col_ind[begin + tid]];
	}
      // if we need to reduce array to get into shared.
      else
	{
	  //printf("in multichunk\n");
	  int i;
	  for (i = begin + tid*chunk ; i < begin + (tid+1)*chunk ; i++)
	    {
	      reduce[tid] += vals[i] * x[col_ind[i]];
	    }
	  
	  // if we have a remainder just use the last thread to get the rest
	  if (tid == (THREADS_PER_BLOCK - 1))
	    {
	      while (i < begin + length)
		{
		  reduce[tid] += vals[i] * x[col_ind[i]];
		  i++;
		}
	    }
	}


      //if (row == 0)
      //	printf("thread %d has %f\n", tid, reduce[tid]);

         
      
      // if reduce is odd number length, and is last thread in block
      if (THREADS_PER_BLOCK % 2 == 1 && tid == THREADS_PER_BLOCK-2)
	reduce[tid] += reduce[tid + 1];
      __syncthreads();

      // do sequential reduction
      for (int s = THREADS_PER_BLOCK/2 ; s > 0 ; s>>=1)
	{
	  if (tid < s)
	    reduce[tid] += reduce[tid + s];
	  __syncthreads();
	}
      
      if (tid == 0)
	b[row] = reduce[tid];
    }
  
  /*
  if (tid == 0 && row == 0)
    {
      printf("\n");
      for (int i = 0 ; i < m ; i++)
	printf("\t%f", b[i]);
      printf("\n");
    }
  */


}



void spmv_gpu(unsigned int* row_ptr, unsigned int* col_ind, double* vals,
              int m, int n, int nnz, double* x, double* b)
{
    // timers
    cudaEvent_t start;
    cudaEvent_t stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float elapsedTime;

    // GPU execution parameters
    // 1 thread block per row
    // 64 threads working on the non-zeros on the same row
    unsigned int blocks = m; 
    //printf("\n\nnumber of blocks: %d\n\n", m);
    unsigned int threads = THREADS_PER_BLOCK; 
    unsigned int shared = threads * sizeof(double);

    dim3 dimGrid(blocks, 1, 1);
    dim3 dimBlock(threads, 1, 1);

    checkCudaErrors(cudaEventRecord(start, 0));
    for(unsigned int i = 0; i < MAX_ITER; i++) {
      cudaDeviceSynchronize();
      spmv_kernel<double><<<dimGrid, dimBlock, shared>>>(row_ptr, col_ind, 
							 vals, m, n, nnz, 
							 x, b);
    }
    checkCudaErrors(cudaEventRecord(stop, 0));
    checkCudaErrors(cudaEventSynchronize(stop));
    checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, stop));
    printf("  Exec time (per itr): %0.8f s\n", (elapsedTime / 1e3 / MAX_ITER));

}
