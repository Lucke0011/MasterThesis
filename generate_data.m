function [data, timelock] = generate_data(signal ,t)
    load('grad.mat')
    n_trials = 1;
    n_samples = length(t);

    data = [];
    data.label = grad.label;
    data.time = cell(1, n_trials);
    for i = 1:n_trials
        data.time{1,i} = t;
    end
    
    % trial
    data.trial = cell(1, n_trials);
    for i = 1:n_trials
        data.trial{1, i} = signal;
    end
    
    data.sampleinfo = [((0:(n_trials-1))*n_samples)+1; ((1:n_trials)*n_samples)]';
    data.trialinfo = ones(n_trials, 1);
    data.grad = grad;

    % for i = 1:length(grad.chantype)
    %     data.grad.chantype{i} = {'MEG'};
    % end

    data = ft_datatype_raw(data);

    timelock        = [];
    timelock.dimord = 'chan_time';
    timelock.avg    = signal;
    timelock.label  = grad.label;
    timelock.time   = t;
    timelock.grad   = grad;

    timelock = ft_datatype_timelock(timelock);
end

