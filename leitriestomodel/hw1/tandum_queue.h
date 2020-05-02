#define Q_LIMIT 100
#define BUSY 1
#define FREE 0


void init_prog(int ss[], int in_q[], float *t, float *t_prev,
	       float next_events[][2], float arrival_mean);

void init_stats(int total_thru[], float total_delay[],
		float in_q_total[], float serv_time_total[]);

void update_t_and_e(float *t, float next_events[][2],
		    int *next_event_type, int *next_event_queue);

void arrive(int neq,float next_events[][2],
	    int ss[], int in_q[], float total_delay[], int total_thru[],
	    float **q_arrive_t, float arrival_mean,
	    float serv_mean[], float t);

void depart(int neq, int in_q[], int ss[],
	    float next_events[][2],
	    float **q_arrive_t,
	    float total_delay[], int total_thru[],
	    float serv_mean[], float t);
  
float expon(float mean);

void update_stats(float t, float *t_prev, int* ss , int *in_q,
		  float *in_q_total, float *serv_time_total);

void printMatrix(float **m, int height, int length);
void print_next_events(float m[][2], int num_s, float T);
