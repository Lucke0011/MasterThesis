%% Data
load('grad.mat')

% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
n_trials = 3;
t = (0:(length(t_trial)-1)) * (1/Fs);

data = [];
%%

% data.label = cell(n_sensors, 1);
% for i = 1:n_sensors
%     data.label{i} = int2str(i);
% end

data.label = grad.label;

%%

data.time = cell(1, n_trials);
for i = 1:n_trials
    data.time{1,i} = t;
end

% trial
data.trial = cell(1, n_trials);
signals = generate_signals(n_trials);
for i = 1:n_trials
    data.trial{1, i} = signals{i};
end

data.sampleinfo = [((0:(n_trials-1))*n_samples)+1; ((1:n_trials)*n_samples)]';
data.trialinfo = ones(n_trials, 1);
data.grad = grad;

%% Test data structure

data_result = ft_datatype_raw(data);

%% Timelocked structure

timelock        = [];
timelock.dimord = 'chan_time';
timelock.avg    = signals{1};
timelock.label  = grad.label;
timelock.time   = t;
timelock.grad   = grad;

%% Test timelock structure

timelock_result = ft_datatype_timelock(timelock);

%% Topoplot
cfg = [];
cfg.grad = grad;
cfg.channel = 'all'; % 'M*' skip the reference channels
cfg.skipscale = 'yes';
cfg.skipcomnt = 'yes';
cfg.projection = 'orthographic';
cfg.width  = 2.0;
cfg.height = 1.5;
layout_orthographic = ft_prepare_layout(cfg);

figure
ft_plot_layout(layout_orthographic)

%%
cfg = [];
cfg.grad = grad;
cfg.channel = 'all';
cfg.skipscale = 'yes';
cfg.skipcomnt = 'yes';
cfg.projection = 'polar';
cfg.width  = 0.2; % the projected sensor positions are dimensionless, this requires some trial and error
cfg.height = 0.15;
layout_polar = ft_prepare_layout(cfg);

% increase the spacing of the channels and shift them a bit
layout_polar.pos(:,1) = layout_polar.pos(:,1) * 1.0;
layout_polar.pos(:,2) = layout_polar.pos(:,2) * 1.2 + 0.02;

figure
ft_plot_layout(layout_polar)

%%

cfg = [];
cfg.layout = layout_polar;

ft_topoplotER(cfg, timelock_result);