using DelimitedFiles
using Distributions

struct queue
    
    number::Int64
    serv_mean::Float64
    # travel time mean
    tt_mean::Float64
    # arrival times in the queue for this server
    q_arrival_t::Array{Float64,1}
    # current travel times to this queue
    travel_time::Array{Float64,1}
    # server status
    free::Bool
    total_thru::Int64
    total_delay::Float64
    # number of people in queue weighted by time intervals
    tw_queue::Float64
    
end

struct event

    time::Float64
    queue_num::Int64
    type::Int64

end
    



expon(mean) = rand(Exponential(mean),1)[1]

function main()

    time_data = readdlm("tandum_queue.in")
    num_queues = Int(time_data[1,1])
    T = time_data[1,2]
    wait_times = time_data[2:end, :]

    # make array of queue structs to store all data about
    # each queue
    queues = Array{queue, 1}(undef, num_queues)
    for i=1:num_queues
        queues[i] = queue(i, wait_times[i,2], wait_times[i,1],
                          [], true, 0, 0.0, 0.0)    
    end

    current_event =  event(expon(queues[1].tt_mean), 1, 1)
    events = []
  
    t = 0.0

    while t < T
        
        if current_event.type == 1
            arrival!(t, queues[current_event.queue], events)
        elseif current_event.type == 2
            departure!(queues, current_event)
        end
        
    end
        
    
end

function arrival!(t, queue, events)

    if queue.number == 1
        push!(events, event(t + expon(queue.tt_mean), 1, 1))
    end

    if queue.free == false
        push!(q_arrival_t, t)
    else
        queue.total_thru += 1
        queue.free = false
        push!(events, event(t + expon(queue.serv_mean), 1, 2))
    end

end




main()
