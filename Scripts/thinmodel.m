function md = thinmodel(md,timesteps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    time = [];
    for i = 1:length(md.results.TransientSolution)
        time = [time md.results.TransientSolution(i).time];
    end

    time_index = [];
    for i = 1:length(timesteps)
        [val,idx]=min(abs(time-timesteps(i))); % Closted timestep to the prescribed dates 
        time_index = [time_index idx];
    end

    time_index = unique(time_index);
    md.results.TransientSolution = md.results.TransientSolution(time_index);

end