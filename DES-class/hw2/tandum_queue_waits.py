from numpy.random import exponential
from numpy.random import uniform


class Queue:
    
    def __init__(self, number, serv_mean, tt_mean):

        self.number = number
        self.serv_mean= serv_mean
        # this is the arrival mean in the case of queue 0, and the uniform range in case of other queues
        self.tt_mean = tt_mean
        # arrival times in queue
        self.q_arrival_t = []
        # current number of people traveling to the queue
        self.num_transit = 0
        self.busy = False
        self.total_thru = 0
        self.total_delay = 0.0
        # number of people in queue weighted by time intervals
        self.tw_queue = 0.0
        self.serv_time = 0.0
        self.total_transit = 0.0

    def __str__(self):

        to_print = "__________________\n" +"Queue Number " + str(self.number) + "\n"

        to_print += "Mean Serving Time: " + str(self.serv_mean) + "\n"
        to_print += "Mean Time to Travel To: " + str(self.tt_mean) + "\n"
        to_print += "Current Number in Queue: " + str(len(self.q_arrival_t)) + "\n"
        to_print += "Current Number Traveling to Queue/Server: " + str(self.num_transit) + "\n"
        to_print += "Servers Current Status: " + str(self.busy) + "\n"
        to_print += "Total Thru Currently: " + str(self.total_thru) + "\n"
        
        return to_print
        
    def update_stats(self, t_period):

        self.tw_queue += len(self.q_arrival_t) * t_period
        self.total_transit += self.num_transit * t_period
        self.serv_time += int(self.busy) * t_period
        
        
class Event:
    
    def __init__(self, t, queue_num, type_e):
        self.t = t
        self.queue_num = queue_num
        self.type_e = type_e  

class Event_List:

    def __init__(self):
        self.event_list = []
        self.total_arrivals = 0

    def __str__(self):

        to_print = "now--"
        for event in self.event_list:
            if event.type_e == 0:
                to_print += "(arrival at " + str(event.queue_num) + " at " + str(round(event.t,2)) +")--"
            if event.type_e == 1:
                to_print += "(departure from " + str(event.queue_num) + " at " + str(round(event.t,2)) +")--"
                                                                                       
        to_print += "(end)"
        return to_print
        
        

    def add_event(self, t, queue_num, type_e):

        if type_e == 0 and queue_num == 0:
            self.total_arrivals += 1
        
        
        self.event_list.append(Event(t, queue_num, type_e))
        self.event_list.sort(key = lambda event: event.t)
        

    def peek_next(self):
        return self.event_list[0]
        
    def pop_next(self):
        return self.event_list.pop(0)




def expon(mean):
    return exponential(scale=mean)

def arrive(t, queue, events):

    # if we are at the first queue
    if queue.number == 0:
        # schedule the next arrival
        #print("inter arrival mean:", queue.tt_mean)
        events.add_event(t + expon(queue.tt_mean), 0, 0)
    else:
        #if this isn't the first queue person is no longer in transit
        queue.num_transit -= 1

    # if they have to wait in line
    if queue.busy == True:
        # add there arrival time
        queue.q_arrival_t.append(t)
        
    # server is free can serv immediatly
    else:
        #print("person immideitly counted")
        # add them to tru
        queue.total_thru += 1
        #print(queue.total_thru)
        # change server status to busy
        queue.busy = True
        # schedule there departure
        events.add_event(t + expon(queue.serv_mean), queue.number, 1)


def departure(t, dep_queue, queues, current_event, events):

    c_num = dep_queue.number
    

    # if no one is waiting after this depature
    if len(dep_queue.q_arrival_t) == 0:
        dep_queue.busy = False

    # if someone is waiting
    else:
        # add their delay
        #print(dep_queue.q_arrival_t[0], t)
        #print(t - dep_queue.q_arrival_t[0])
        dep_queue.total_delay += t - dep_queue.q_arrival_t.pop(0)
        # add them to thru
        #print("person waiting counted")
        dep_queue.total_thru += 1
        #print(dep_queue.total_thru)
        # schedule there departure
        #print("serv mean: ", dep_queue.serv_time)
        events.add_event(t + expon(dep_queue.serv_mean), c_num, 1)

    # if this is not the last queue
    if dep_queue.number != len(queues) - 1:
        # have this customer arrive in next queue after travel time
        low = queues[c_num+1].tt_mean[0]
        high = queues[c_num+1].tt_mean[1]
        events.add_event(t + uniform(low=low, high=high), c_num + 1, 0)
        queues[c_num+1].num_transit += 1
        
def write_stats(queues, fout, t):

    
    fout.write("\n\nAverage delay in queues:\n\n")
    for i in range(len(queues)):
        fout.write("Queue %d\t\t" % i)
    fout.write("\n\n")
    for queue in queues:
        #print(queue.total_delay, queue.total_thru)
        avg_delay = queue.total_delay/queue.total_thru
        fout.write(str(avg_delay) + "\t\t")
    
    fout.write("\n\n")
    fout.write("\n\nAverage number in queues:\n\n")
    for i in range(len(queues)):
        fout.write("Queue %d\t\t" % i)
    fout.write("\n\n")
    for queue in queues:
        avg_num_in_q = queue.tw_queue/t
        fout.write(str(avg_num_in_q) + "\t\t")
    fout.write("\n\n")
    fout.write("\n\nServer utilization:\n\n")
    for i in range(len(queues)):
        fout.write("Queue %d\t\t" % i)
    fout.write("\n\n")
    for queue in queues:
        serv_util = queue.serv_time/t
        fout.write(str(serv_util) + "\t\t")
    fout.write("\n\n")
    fout.write("\n\nAverage Number in Transit to :\n\n")
    for i in range(len(queues)):
        fout.write("Queue %d\t\t" % i)
    fout.write("\n\n")
    for queue in queues:
        avg_transit = queue.total_transit/t
        fout.write(str(avg_transit) + "\t\t")
    fout.write("\n\n")


    
        
def main():
    
    fin = open("tandum_queue.in", "r")
    f_line = [int(i) for i in fin.readline().strip("\n").split('\t')]
    num_queues = f_line[0]
    T = f_line[1]
    tt_uni = [float(i) for i in fin.readline().strip("\n").split('\t')]

    serv_means = [float(i) for i in fin.read().split('\n')]
    
    
    fout = open("tandum_queue.out", "w")

    fout.write("System of " + str(num_queues) + " servers in tandem running for " +  str(T) + " minutes\n\n")
    fout.write("Mean arrival time: " + str(serv_means[0]) + '\n')
    fout.write("Minimum and Maximum transit times: " + str(tt_uni) + '\n')
    for i in range(1, num_queues+1):
        fout.write("Mean service time of queue " + str(i-1) + ": " + str(serv_means[i]) + '\n') 


    for rep in range(10):

        fout.write("\nRepatition number " + str(rep))
        # making the data for each queue
        queues = []
        for i in range(num_queues):
            if i == 0:
                q = Queue(i, serv_means[1] , serv_means[0])
            else:    
                q = Queue(i, serv_means[i+1] , tt_uni)
            queues.append(q)

        t = 0.0
        t_prev = 0.0  
        # making the events list and scheduling first event

        events = Event_List()
        events.add_event(t + exponential(scale=queues[0].tt_mean), 0, 0)
        
        while t < T:

            current_event = events.pop_next()
            t = current_event.t
            t_inter = t - t_prev
        
            # update ongoing stats
            for queue in queues:
                queue.update_stats(t_inter)
            
            if current_event.type_e == 0:
                arrive(t, queues[current_event.queue_num], events)
            
            if current_event.type_e == 1:
                departure(t, queues[current_event.queue_num], queues, current_event, events)
            
            t_prev = t

        write_stats(queues, fout, t)
    
            
main()
