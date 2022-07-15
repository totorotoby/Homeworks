#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "common.h"


void usage(int argc, char** argv);
void verify(int* sol, int* ans, int n);
void prefix_sum(int* src, int* prefix, int n);
void prefix_sum_p1(int* src, int* prefix, int n);
void prefix_sum_p2(int* src, int* prefix, int n);


int main(int argc, char** argv)
{
    // get inputs
    uint32_t n = 1048576;
    unsigned int seed = time(NULL);
    if(argc > 2) {
        n = atoi(argv[1]); 
        seed = atoi(argv[2]);
    } else {
        usage(argc, argv);
        printf("using %"PRIu32" elements and time as seed\n", n);
    }



    // set up timers
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();

    // files to write to
    FILE *fout1 = fopen("tests", "a");
    


    // set up data 
    int* prefix_array = (int*) AlignedMalloc(sizeof(int) * n);  
    int* input_array = (int*) AlignedMalloc(sizeof(int) * n);
    srand(seed);
    
	  
    for(int j = 0; j < n; j++)
      input_array[j] = rand() % 100;

       
    // execute serial prefix sum and use it as ground truth
    start_t = ReadTSC();
    prefix_sum(input_array, prefix_array, n);
    end_t = ReadTSC();
    printf("Time to do O(N-1) prefix sum on a %"PRIu32" elements: %g (s)\n", 
	   n, ElapsedTime(end_t - start_t));

    fprintf(fout1, "%g,", ElapsedTime(end_t - start_t));
		
    int* input_array1 = (int*) AlignedMalloc(sizeof(int) * n);  
    int* prefix_array1 = (int*) AlignedMalloc(sizeof(int) * n);
	
    memcpy(input_array1, input_array, sizeof(int) * n);
    
    start_t = ReadTSC();
    prefix_sum_p1(input_array1, prefix_array1, n);
    end_t = ReadTSC();
    printf("Time to do O(NlogN) //prefix sum on a %"PRIu32" elements: %g (s)\n",
			 n, ElapsedTime(end_t - start_t));
    verify(prefix_array, prefix_array1, n);
	
    fprintf(fout1, "%g,", ElapsedTime(end_t - start_t));
	


    // execute parallel prefix sum which uses a 2(N-1) algorithm
    memcpy(input_array1, input_array, sizeof(int) * n);
    memset(prefix_array1, 0, sizeof(int) * n);
    start_t = ReadTSC();
    prefix_sum_p2(input_array1, prefix_array1, n);
    end_t = ReadTSC();
    printf("Time to do 2(N-1) //prefix sum on a %"PRIu32" elements: %g (s)\n", 
	   n, ElapsedTime(end_t - start_t));
    verify(prefix_array, prefix_array1, n);


    fprintf(fout1, "%g\n", ElapsedTime(end_t - start_t));

	
    // free memory
    AlignedFree(prefix_array);
    AlignedFree(input_array);
    AlignedFree(input_array1);
    AlignedFree(prefix_array1);

      
    
    return 0;
}

void usage(int argc, char** argv)
{
    fprintf(stderr, "usage: %s <# elements> <rand seed>\n", argv[0]);
}


void verify(int* sol, int* ans, int n)
{
    int err = 0;
    for(int i = 0; i < n; i++) {
      //printf("%d, %d\n", sol[i], ans[i]);
        if(sol[i] != ans[i]) {
            err++;
        }
    }
    
    if(err != 0) {
        fprintf(stderr, "There was an error: %d\n", err);
    } else {
        fprintf(stdout, "Pass\n");
    }
}

void prefix_sum(int* src, int* prefix, int n)
{
    prefix[0] = src[0];
    for(int i = 1; i < n; i++) {
        prefix[i] = src[i] + prefix[i - 1];
    }
}

void prefix_sum_p1(int* src, int* prefix, int n)
{

  int i;
  int j;
  // number of levels in tree
  int tree_iter = log2(n);
  // stride between elements to be added
  int stride = 1;

  // loop over number of passes on array combining 2 elements
  for (i = 0 ; i < tree_iter ; i++)
    {
      // loop over elements in array
      #pragma omp parallel for private(j)
      for (j = 0 ; j < n; j++)
	{
	  // update indices that are not already done on current iteration
	  if (j >= stride)
	    prefix[j] = src[j] + src[j-stride];
	  
	  // otherwise keep the done values around to be used later
	  else
	    prefix[j] = src[j];
	}
      
      // in order to avoid WARs need to copy back to orginal array to read from
      #pragma omp parallel for private(j)
      for (j = 0 ; j < n; j++)
	src[j] = prefix[j];
      
      // update the stride which is 2^i
      stride = stride * 2;
    }
}


void prefix_sum_p2(int* src, int* prefix, int n)
{
  
  int i;
  int j;

  // number of levels in tree
  int tree_iter = log2(n);
  // stride to read back to
  int rstride = 1;
  // stride between elements on current level
  int cstride = 2;

  //copy src into prefix
  #pragma omp parallel for private(j)
  for (j = 0 ; j < n ; j++)
    prefix[j] = src[j];

  
  // reduce phase
  for (i = 0 ; i < tree_iter ; i++)
    {
      //traverse nodes a current level to update to next level
      #pragma omp parallel for private(j)
      for (j = 0 ; j < n ; j = j + cstride)
	// traverse tree upward by adding left node and right node togther
	// storing result in parent node stored at the index of previous right node.
	prefix[j + cstride - 1] = prefix[j + cstride - 1] + prefix[j + rstride - 1];
      rstride = 2 * rstride;
      cstride = 2 * cstride;
    }


  // stride between elements on current level to be read
  cstride = n/2;
  // stride between elements on next level to write to
  int wstride = n/4;
  // down add phase
  for (i = 0 ; i < tree_iter-1 ; i++)
    {
      //printf("on level %d:\n", i);
      #pragma omp parallel for private(j)
      for (j = cstride ; j < n ; j = j + cstride)
	prefix[j + wstride - 1] += prefix[j-1];
	  
      cstride = cstride/2;
      wstride = wstride/2;   
    }
}
