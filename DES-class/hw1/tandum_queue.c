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
  fout = fopen("tandum_queue_problem_3-1.out", "w");
  fscanf(fin, "%d %f\n", &num_s, &T);

  float serv_mean[num_s];
  
  fscanf(fin, "%f %f", &arrival_mean, &serv_mean[0]);
  for (i = 0 ; i < num_s - 1 ; i++)
    fscanf(fin, "\n%f", &serv_mean[i + 1]);
    
  fprintf(fout, "System of %d servers in tandem running %.2f minutes\n\n",
	  num_s, T);

  
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

  // repatitions loop
  for (int j = 0 ; j < 10 ; j++)
    {
      
      init_prog(ss, in_q, &t, &t_prev, next_events, arrival_mean);

  
      init_stats(total_thru, total_delay,
		 in_q_total,
		 serv_time_total);
      
      fprintf(fout, "\n\n Repatition #%d\n\n", j+1);
      fprintf(fout, "Mean arrival time: %.2f\n\n", arrival_mean);
      for (i = 0 ; i < num_s ; i++)
	fprintf(fout, "Mean service time of server %d: %.2f\n\n",
		i+1, serv_mean[i]);
      
  
      
      /*start simulating and stop when end time is reached*/
      while ( t < T )
	{

      
	  update_t_and_e(&t, next_events, &next_event_type, &neq);

	  update_stats(neq, t, &t_prev, ss,
		       in_q,
		       in_q_total,
		       serv_time_total);

      
	  /*print whats happening */
	  //************************************************************
      
	  printf("\n*********************************************************\n");
	  printf("Current Time: %f\n", t);
	  printf("Current Action:\n");
	  if (next_event_type == 0){
	    printf("customer is about to arrive at queue %d\n\n", neq);
	  }
	  if (next_event_type == 1){
	    printf("customer is about to depart from queue %d\n\n", neq);
	  }
	  
	  printf("Next Events Scheduled:\n\n");
	  print_next_events(next_events, num_s);
	  printf("\n");
	  printf("Number of people currently in queues:\n");
	  print_queue_stat_i(in_q, num_s);
	  printf("Number of people total through queues:\n");
	  print_queue_stat_i(total_thru, num_s);
	  printf("The Total delay in each queue:\n");
	  print_queue_stat_f(total_delay, num_s);
	  printf("The Total time server is busy:\n");
	  print_queue_stat_f(serv_time_total, num_s);
	  printf("The Total in queue integrate by time:\n");
	  print_queue_stat_f(in_q_total, num_s);
      
      //**************************************************************************


      
	  if (next_event_type == 0){
	    arrive(neq, next_events, ss, in_q,
		   total_delay, total_thru,
		   q_arrive_t, arrival_mean, serv_mean, t, fout);
	  }

	  if (next_event_type == 1){
	
	    depart(neq, in_q, ss,
		   next_events, q_arrive_t,
		   total_delay, total_thru,
		   serv_mean, t);
	  }
	}

      write_stats(fout, num_s, t, total_delay, total_thru,
		  in_q_total, serv_time_total);

      //problem three decreasing
      //arrival_mean =arrival_mean - .0575;
      //serv_mean[0] = serv_mean[0] - .0275;
      
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

void update_stats(int neq, float t, float *t_prev, int ss[] , int in_q[],
		  float in_q_total[], float serv_time_total[])
{
  float time_interval;

  time_interval = t - (*t_prev);
  /* update for the next time step */
  *t_prev = t;

  serv_time_total[neq] += ss[neq] * time_interval;
  in_q_total[neq] += in_q[neq] * time_interval;
 
}

void arrive(int neq,float next_events[][2],
	    int ss[], int in_q[], float total_delay[],
	    int total_thru[],float **q_arrive_t,
	    float arrival_mean, float serv_mean[], float t, FILE *fout)
{


  /* If this is the first queue we need to schedule the next arrive to the system */
  if(neq == 0)
    next_events[neq][0] = t + expon(arrival_mean);
  /* If this isnt the first queue then arrrivals can only be triggered by departures
     so make arrival impossible untill next departure */
  if (neq >= 1)
    next_events[neq][0] = FLT_MAX;
  

  /* if someone is already being served need to join the queue */
  if (ss[neq] == BUSY)
    {
      in_q[neq]++;
      if (in_q[neq] > Q_LIMIT)
	{
	  fprintf(fout, "Overflow in queue 0 at time %f\n", t);
	  exit(2);
	}
      
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
      /* compute the delay of the new person */
      delay = t - q_arrive_t[neq][in_q[neq]];
      /*printf("Time is: %f    arrival time is: %f    delay is: %f\n",
	t, q_arrive_t[neq][in_q[neq]], delay);*/
      total_delay[neq] += delay;
     
      /* There is one less in the queue now because we took a new person */
      in_q[neq]--;
      /* increment total peopel thru */
      total_thru[neq]++;

      /* add departure time for this new person */
      next_events[neq][1] = t + expon(serv_mean[neq]);

      
      /* queue needs to be moved up */
      for (i = 1 ; i <= in_q[neq] ; i++)
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

void write_stats(FILE *fout, int num_s, float t,
		 float total_delay[], int total_thru[], float in_q_total[],
		 float serv_time_total[])
{

  int i;
  
  fprintf(fout, "\n\nAverage delay in queues:\n\n");
  for (i = 0 ; i < num_s ; i++)
    fprintf(fout, "Queue %d\t\t", i);
  fprintf(fout, "\n\n");
  //for (i = 0 ; i < num_s*15 ;i++)
    //fprintf(fout, "-");
  fprintf(fout, "\n\n");
  for (i = 0 ; i < num_s ; i++)
    fprintf(fout, "%2.5f\t\t" , total_delay[i]/total_thru[i]);
  fprintf(fout, "\n\n");
     fprintf(fout, "\n\nAverage number in queues:\n\n");
  for (i = 0 ; i < num_s ; i++)
     fprintf(fout, "Queue %d\t\t", i);
   fprintf(fout, "\n\n");
   //for (i = 0 ; i < num_s*15 ;i++)
    // fprintf(fout, "-");
   fprintf(fout, "\n\n");
  for (i = 0 ; i < num_s ; i++)
     fprintf(fout, "%2.5f\t\t" , in_q_total[i]/t);
   fprintf(fout, "\n\n");
    fprintf(fout, "\n\nServer utilization:\n\n");
  for (i = 0 ; i < num_s ; i++)
     fprintf(fout, "Queue %d\t\t", i);
   fprintf(fout, "\n\n");
   //for (i = 0 ; i < num_s*15 ;i++)
    //     fprintf(fout, "-");
   fprintf(fout, "\n\n");
  for (i = 0 ; i < num_s ; i++)
     fprintf(fout, "%2.5f\t\t" , serv_time_total[i]/t);
   fprintf(fout, "\n\n");
    

}




void printMatrix(float **m, int height, int length){

  for (int i = 0 ; i < height ; i++){
    //printf("\n");
    for (int j = 0 ; j < length ; j++){
      printf("%2.2e,", m[i][j]);
    }
    printf("\n");
  }
}

void print_next_events(float m[][2], int num_s){

  for (int i = 0 ; i < num_s ; i++){
    for (int j = 0 ; j < 2; j++){
      if (j == 0 && m[i][j] != FLT_MAX)
	printf("Next arrival at queue %d is at %2.2f\n", i, m[i][j]);
      if (j == 0 && m[i][j] == FLT_MAX)
	printf("Next arrival at queue %d is not scheduled\n", i);
      if (j == 1 && m[i][j] != FLT_MAX)
	printf("Next departure at queue %d is at %2.2f\n", i, m[i][j]);
      if (j == 1 && m[i][j] == FLT_MAX)
	printf("Next departure at queue %d is not scheduled\n", i);
	
    }
    printf("\n");
  }
}

void print_queue_stat_i(int stat[], int num_s)
{

  int i;
  
  for (i = 0 ; i < num_s ; i++)
    printf("Queue %d\t\t", i);
  printf("\n");
  for (i = 0 ; i < num_s*15 ;i++)
    printf("-");
  printf("\n");
  for (i = 0 ; i < num_s ; i++)
    printf("%d\t\t" , stat[i]);
  printf("\n");

}

void print_queue_stat_f(float stat[], int num_s)
{

  int i;
  
  for (i = 0 ; i < num_s ; i++)
    printf("Queue %d\t\t", i);
  printf("\n");
  for (i = 0 ; i < num_s*15 ;i++)
    printf("-");
  printf("\n");
  for (i = 0 ; i < num_s ; i++)
    printf("%2.2f\t\t" , stat[i]);
  printf("\n");

}




float expon(float mean){return -mean * log(lcgrand(1));}
