#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lcgrand.h"
#include <float.h>
#include "tandum_queue.h"

int num_s;

int main()
{
  FILE *fin, *fout;

  /*Model input parameters */

  float T, arrival_mean;
  int i;
  
  fin = fopen("tandum_queue.in", "r");
  fout = fopen("tandum_queue.out", "w");
  fscanf(fin, "%d %f\n", &num_s, &T);

  float serv_mean[num_s];
  
  fscanf(fin, "%f %f", &arrival_mean, &serv_mean[0]);
  for (i = 0 ; i < num_s - 1 ; i++)
    fscanf(fin, "\n%f", &serv_mean[i + 1]);
    
  fprintf(fout, "System of %d servers in tandem running %.2f minutes\n\n",
	  num_s, T);
  fprintf(fout, "Mean arrival time: %.2f\n", arrival_mean);
  for (i = 0 ; i < num_s ; i++)
    fprintf(fout, "Mean service time of server %d: %.2f\n",
	    i+1, serv_mean[i]);
  
  /* program data */
  
  int ss[num_s],
    in_q[num_s],
    next_event_type,
    /* neq stands for next event queue */
    neq;
  
  float
    t,
    t_prev,
    next_events[num_s][2];

  /* making arrive times on the heap so that queue length can be big */
  float **q_arrive_t = malloc(sizeof(float*) * num_s);
  for (i = 0 ; i < num_s ; i++)
    {
      q_arrive_t[i] = malloc(sizeof(float) * (Q_LIMIT + 1));
    }
  
  /* statistical data */
  int total_thru[num_s];

  float total_delay[num_s],
    in_q_total[num_s],
    serv_time_total[num_s];
  
  
  init_prog(ss, in_q, &t, &t_prev, next_events, arrival_mean);

  
  init_stats(total_thru, total_delay,
	     in_q_total,
	     serv_time_total);
  
  
  /*start simulating and stop when end time is reached*/
  while ( t < T )
    {

      
      update_t_and_e(&t, next_events, &next_event_type, &neq);

      update_stats(t, &t_prev, &ss[neq],
		   &in_q[neq],
		   &in_q_total[neq],
		   &serv_time_total[neq]);

      /*print whats happening */
      printf("\n------------------------------------\n");
      printf("Current Time: %f\n", t);
      if (next_event_type == 0){
	printf("customer is about to arrive at queue %d\n", neq);
      }
      if (next_event_type == 1){
	printf("customer is about to depart from queue %d\n", neq);
      }

       
      if (next_event_type == 0){
	  arrive(neq, next_events, ss, in_q,
		 total_delay, total_thru,
		 q_arrive_t, arrival_mean, serv_mean, t);
      }

      if (next_event_type == 1){
	printf("in queue array: %d %d\n", in_q[0], in_q[1]);
      depart(neq, in_q, ss,
	     next_events, q_arrive_t,
	     total_delay, total_thru,
	     serv_mean, t);
      }
    }
  
  fclose(fin);
  fclose(fout);
  free(q_arrive_t);
}  

void init_prog(int ss[], int in_q[], float *t, float *t_prev,
	       float next_events[][2], float arrival_mean)
{

  int i;
  for(i = 0 ; i < num_s ; i++)
    {
      /* all servers start free */
      ss[i] = FREE;
      /* no one is in any queues yet */
      in_q[i] = 0;
      /* force any departure to be impossible on first step */
      next_events[i][1] = FLT_MAX;
      /* force any arrivals to be impossible not at the first server */
      if (i != 0)
	next_events[i][0] = FLT_MAX;
    }
  
  *t = 0.0;
  next_events[0][0] = *t + expon(arrival_mean);
}

void init_stats(int total_thru[], float total_delay[],
		float in_q_total[],
		float serv_time_total[])
{

  int i;
  for (i = 0 ; i < num_s ; i++)
    {
      total_thru[i] = 0;
      total_delay[i] = 0.0;
      in_q_total[i] = 0.0;
      serv_time_total[i] = 0.0;
    }
}

void update_t_and_e(float *t, float next_events[][2],
		    int *next_event_type, int *next_event_queue)
{

  int i;
  int j;
  float minimum = FLT_MAX;

  printMatrix_2(next_events, num_s);

  
  for (i = 0 ; i < num_s ; i ++)
    {
      for (j = 0 ; j < 2 ; j ++)
	{
	  if (minimum > next_events[i][j])
	    {
	      minimum = next_events[i][j];
	      *next_event_queue = i;
	      /* 0 is an arrival in the queue, 1 is a departure in the queue */
		*next_event_type = j;
	    }
	}
    }
  *t = minimum;
  
}

void update_stats(float t, float *t_prev, int* ss , int *in_q,
		  float *in_q_total, float *serv_time_total)
{
  int time_interval;

  time_interval = t - (*t_prev);
  /* update for the next time step */
  *t_prev = t;

  
  *serv_time_total += (*ss) * time_interval;
  *in_q_total += (*in_q) * time_interval;
 
}

void arrive(int neq,float next_events[][2],
	    int ss[], int in_q[], float total_delay[],
	    int total_thru[],float **q_arrive_t,
	    float arrival_mean, float serv_mean[], float t)
{


  /* If this is the first queue we need to schedule the next arrive to the system */
  if(neq == 0)
    next_events[neq][0] = t + expon(arrival_mean);
  /* If this isnt the first queue then arrrivals can only be triggered by departures
     so make arrival impossible untill next departure */
  if (neq >= 1)
    next_events[neq][0] = FLT_MAX;
  

  /* if someone is already being served need to join the queue */
  printf("sever status: %d\n", ss[neq]);
  if (ss[neq] == BUSY)
    {
      in_q[neq]++;
      if (in_q[neq] > Q_LIMIT)
	{
	  printf("Overflow in queue 0 at time %f\n", t);
	  exit(2);
	}
      printf("number currently in queue: %d\n", in_q[neq]);
      q_arrive_t[neq][in_q[neq]] = t;
    }

  /* if no one is being served then we can just serve the person */
  else{
    
    total_delay[neq] += 0;
    total_thru[neq] += 1;

    ss[neq] = BUSY;

    next_events[neq][1] = t + expon(serv_mean[neq]);
  }
}

void depart(int neq, int in_q[], int ss[],
	    float next_events[][2],
	    float **q_arrive_t,
	    float total_delay[], int total_thru[],
	    float serv_mean[], float t)
{

  float delay;
  int i;

  /* case where no one is in the queue and the last person departs */
  if (in_q[neq] == 0)
    {
      ss[neq] = FREE;
      next_events[neq][1] = FLT_MAX;
    }

  /* case where someone is departing and we need to schedule a new person in */
  else
    {
      /* There is one less in the queue now because we took a new person */
      fprintf(stderr, "number of people in queue: %d\n", in_q[neq]);
      in_q[neq]--;
      /* compute the delay of the new person */
      fprintf(stderr, "number of people in queue: %d\n", in_q[neq]);
      fprintf(stderr, "queue number: %d\n", neq);
      delay = t - q_arrive_t[neq][in_q[neq]];
      total_delay[neq] += delay;
      /* increment total peopel thru */
      total_thru[neq]++;

      /* add departure time for this new person */
      next_events[neq][1] = t + expon(serv_mean[neq]);

      
      /* queue needs to be moved up */
      for (i = 1 ; i <= in_q[neq] ; i++)
	printf("Howdy!\n");
	q_arrive_t[neq][i] = q_arrive_t[neq][i+1];
    }

  /* if this is not the last queue then we need to move the customer to the next
     queue in the  line */
  if (neq != (num_s - 1))
    {
      /* we can do this by just scheduling a new arrival at the next queue for
	 right now */    
      next_events[neq+1][0] = t;
    }
  
}


void printMatrix(float **m, int height, int length){

  for (int i = 0 ; i < height ; i++){
    //printf("\n");
    for (int j = 0 ; j < length ; j++){
      printf("%.2f,", m[i][j]);
    }
    printf("\n");
  }
}

void printMatrix_2(float m[][2], int height){

  for (int i = 0 ; i < height ; i++){
    for (int j = 0 ; j < 2; j++){
      printf("%.2f,", m[i][j]);
    }
    printf("\n");
  }
}

float expon(float mean){return -mean * log(lcgrand(1));}
