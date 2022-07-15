#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "common.h"
#include <inttypes.h>
#include <time.h>


#define PI 3.1415926535


void usage(int argc, char** argv);
double calcPi_Serial(int num_steps);
double calcPi_P1(int num_steps);
double calcPi_P1_threads(int num_steps, int nthreads);
double calcPi_P2(int num_steps);
double calcPi_monte_serial(int num_steps, double seed);
double calcPi_monte_1(int num_steps, double seed);
double calcPi_monte_2(int num_steps, double seed);


int main(int argc, char** argv)
{

  FILE *fout1 = fopen("graph_data", "w");
  FILE *fout2 = fopen("convergence", "w");
  FILE *fout3 = fopen("num_threads", "w");


  fprintf(stdout, "The first 10 digits of Pi are %0.10f\n", PI);
    // set up timer
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();
  
    //seed for random number generator
    double seed = time(0);

    // loop over number of samples for each method
    for (int i = 0 ; i < 7 ; i++){

      // number of samples as powers of 10
      uint32_t num_steps = 10 * pow(10,i);

      if(argc > 1) {
	num_steps = atoi(argv[1]);
      } else {
	usage(argc, argv);
	printf("using %"PRIu32"\n", num_steps);
      }

      // calculate in serial 
      start_t = ReadTSC();
      double Pi0 = calcPi_Serial(num_steps);
      end_t = ReadTSC();
      printf("Time to calculate Pi serially with %"PRIu32" steps is: %g\n",
	     num_steps, ElapsedTime(end_t - start_t));
      printf("Pi is %0.10f\n", Pi0);

      fprintf(fout1, "%.10e,", ElapsedTime(end_t - start_t));
      fprintf(fout2, "%0.10lf & ", Pi0);
    
      // calculate in parallel with reduce
      calcPi_P1(num_steps);
      start_t = ReadTSC();
      double Pi1 = calcPi_P1(num_steps);
      end_t = ReadTSC();
      printf("Time to calculate Pi in // with %"PRIu32" steps is: %g\n",
	     num_steps, ElapsedTime(end_t - start_t));
      printf("Pi is %0.10f\n", Pi1);

      fprintf(fout1, "%.10e,", ElapsedTime(end_t - start_t));
      fprintf(fout2, "%0.10lf & ", Pi1);
      
      // calculate in parallel with atomic add
      calcPi_P2(num_steps);
      start_t = ReadTSC();
      double Pi2 = calcPi_P2(num_steps);
      end_t = ReadTSC();
      printf("Time to calculate Pi in // + atomic with %"PRIu32" steps is: %g\n",
	     num_steps, ElapsedTime(end_t - start_t));
      printf("Pi is %0.10f\n", Pi2);

      fprintf(fout1, "%.10e,", ElapsedTime(end_t - start_t));
      fprintf(fout2, "%0.10lf & ", Pi2);
      
      // calculate monte-carlo serial
      start_t = ReadTSC();
      double Pi_monte = calcPi_monte_serial(num_steps, seed);
      end_t = ReadTSC();
      printf("Time to calculate Pi in serial monte with %"PRIu32" steps is: %g\n",
	     num_steps, ElapsedTime(end_t - start_t));
      printf("Pi is %0.10f\n", Pi_monte);

      fprintf(fout1, "%.10e,", ElapsedTime(end_t - start_t));
      fprintf(fout2, "%0.10lf & ", Pi_monte);
      
      // calculate monte-carlo in parallel with reduce
      start_t = ReadTSC();
      calcPi_monte_1(num_steps, seed);
      double Pi_monte_1 = calcPi_monte_1(num_steps, seed);
      end_t = ReadTSC();
      printf("Time to calculate Pi in parralell for with reduction monte with %"PRIu32" steps is: %g\n",
	     num_steps, ElapsedTime(end_t - start_t));
      printf("Pi is %0.10f\n", Pi_monte_1);

      fprintf(fout1, "%.10e,", ElapsedTime(end_t - start_t));
      fprintf(fout2, "%0.10lf & ", Pi_monte_1);
      
      // calculate monte-carlo in parallel with atomic add
      start_t = ReadTSC();
      calcPi_monte_2(num_steps, seed);
      double Pi_monte_2 = calcPi_monte_2(num_steps, seed);
      end_t = ReadTSC();
      printf("Time to calculate Pi in parralell for with atomic monte with %"PRIu32" steps is: %g\n",
	     num_steps, ElapsedTime(end_t - start_t));
      printf("Pi is %0.10f\n", Pi_monte_2);

      fprintf(fout1, "%.10e\n", ElapsedTime(end_t - start_t));
      fprintf(fout2, "%0.10lf \\\\\n\\hline\n", Pi_monte_2);
    }

    
    // test quad parrallel with different numbers of threads
    for (int i=1 ; i <= 64 ; i++){
      start_t = ReadTSC();
      calcPi_P1_threads(1000000, i);
      end_t = ReadTSC();
      fprintf(fout3, "%.10e\n", ElapsedTime(end_t - start_t));
    }
  
  fclose(fout1);
  fclose(fout2);
  fclose(fout3);
   
  return 0;
}

void usage(int argc, char** argv)
{
    fprintf(stdout, "usage: %s <# steps>\n", argv[0]);
}

double calcPi_Serial(int num_steps)
{

  // distance between points to be evaluated (width)
  double delta_x = 2.0/num_steps;
  double pi = 0.0;
  double x;
  double h;
  
  for (int i = 0 ; i < num_steps ; i++){
    // current x coordinate in [-1,1]
    x = -1.0 + i*delta_x;
    // height of current "bar"
    h = sqrt(1 - x*x);
    // add area of current bar to pi estimate
    pi += delta_x * h;
  }

  
  
  return 2*pi;
}

double calcPi_P1(int num_steps)
{

  // distance between points to be evaluated (width)
  double delta_x = 2.0/num_steps;
  double pi = 0.0;
  
#pragma omp parallel for reduction(+:pi)
  for(int i = 0; i < num_steps; i++) {
    // current x coordinate in [-1,1]
    double x = -1.0 + i*delta_x;
    // height of current "bar"
    double h = sqrt(1 - x*x);
    // add area of current bar to pi estimate
    pi += delta_x * h;
  }
  
    return 2*pi;
}


// same as calcPi_P1 but now number of threads is a parameter
double calcPi_P1_threads(int num_steps, int nthreads)
{
  // distance between points to be evaluated (width)
  double delta_x = 2.0/num_steps;
  double pi = 0.0;
  
#pragma omp parallel for num_threads(nthreads) reduction(+:pi)
  for(int i = 0; i < num_steps; i++) {
    // current x coordinate in [-1,1]
    double x = -1.0 + i*delta_x;
    // height of current "bar"
    double h = sqrt(1 - x*x);
    // add area of current bar to pi estimate
    pi += delta_x * h;
  }
  
    return 2*pi;
}


double calcPi_P2(int num_steps)
{

  // distance between points to be evaluated (width)
  double delta_x = 2.0/num_steps;
  double pi = 0.0;
  
#pragma omp parallel for 
  for (int i = 0; i < num_steps; i++) {
    // current x coordinate in [-1,1]
    double x = -1.0 + i*delta_x;
    // height of current "bar"
    double h = sqrt(1 - x*x);
    // add area of current bar to pi estimate
    #pragma omp atomic
    pi += delta_x * h;
  }
  
    return 2*pi;
}


double calcPi_monte_serial(int num_steps, double seed){

  srand(1);
  
  double x;
  double y;
  
  int inside = 0;

  for (int i = 0 ; i < num_steps ; i++){

    // x and y coorindates of monte sample are independent uniforms between -1 and 1
    x = (double)rand()/RAND_MAX*2.0 - 1.0;
    y = (double)rand()/RAND_MAX*2.0 - 1.0;
    
    //printf("x coord: %f\n", x);
    //printf("y coord: %f\n", y);
    if (sqrt(x*x + y*y) <= 1.0){
      //printf("hello\n");
      inside ++;
    }
  }

  return 4.0*inside/num_steps;
  
}


double calcPi_monte_1(int num_steps, double seed){

  srand(seed);
  int inside = 0;

#pragma omp parallel for reduction(+:inside)
  for (int i = 0 ; i < num_steps ; i++){

    // x and y coorindates of monte sample are independent uniforms between -1 and 1
    double x = (double)rand()/RAND_MAX*2.0 - 1.0;
    double y = (double)rand()/RAND_MAX*2.0 - 1.0;
    
    if (sqrt(x*x + y*y) <= 1.0){
      inside ++;
    }
  }

  return 4.0*inside/num_steps;
  
}


double calcPi_monte_2(int num_steps, double seed){

  srand(seed);
  int inside = 0;

#pragma omp parallel for
  for (int i = 0 ; i < num_steps ; i++){

    // x and y coorindates of monte sample are independent uniforms between -1 and 1
    double x = (double)rand()/RAND_MAX*2.0 - 1.0;
    double y = (double)rand()/RAND_MAX*2.0 - 1.0;

    
    //printf("par x coord: %f\n", x);
    //printf("par y coord: %f\n", y);
    if (sqrt(x*x + y*y) <= 1.0){
      //printf("hello\n");
      #pragma omp atomic
      inside++;
    }
  }

  //printf("inside: %d", inside);
  return 4.0*inside/num_steps;
  
}



