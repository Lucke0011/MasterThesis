

load('sourcemodelT.mat')
load('headmodels.mat')
load('grad.mat')

sourcemodelT.vnn = normals(sourcemodelT.pos,sourcemodelT.tri,'vertex');

% Parameters
L = 10000;              % Length of signal (number of samples)
nr_sensors = 123;
Fs = 1000; % Hz
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
n_trials = 1;
t = (0:(n_trials*length(t_trial)-1)) * (1/Fs);
n_samples = length(t);
i_signal = 9304; % source in 


%% Leadfield for brain
cfg.grid           = sourcemodelT;
cfg.headmodel      = headmodels.headmodel_meg;
cfg.grad           = grad;
leadfield_prepared = ft_prepare_leadfield(cfg); % compute leadfield
n_sources = length(leadfield_prepared.leadfield);
n_sensors = length(leadfield_prepared.label);

%% Brain signal

lf_opm = zeros(n_sensors,n_sources);
for i = 1:n_sources
    vnnT = sourcemodelT.vnn(i,:)'; % transposed from 1x3 to 3x1 to match the leadfield 123x3
    lf_opm(:,i) = leadfield_prepared.leadfield{i} * vnnT; % 123x1
end

% B = LF * q
B = [];

t_trial = -0.5:(1/Fs):(0.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

q_signal = brain_signal(t);
q = zeros(n_sources, 1000);
q(i_signal, :) = q_signal;
B.brain_signal = lf_opm*q; % 123x1000

figure;

subplot(3,1,1);
plot(t, q_signal)
title('Brain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, B.brain_signal)
title('Magnetic field of Brain Signal');
xlabel('Time (s)');
ylabel('T (Tesla)');
grid on;

t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

q_signal = brain_signal(t);
q = zeros(n_sources, n_samples);
q(i_signal, :) = q_signal;
B.brain_signal = lf_opm*q; % 123x1000

subplot(3,1,3);
[pxx_signal, f] = pwelch(B.brain_signal', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of Brain Signal');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;


%% Sensor noise for each sensor 123x1000

B.sensor_noise = zeros(nr_sensors, L);
for i = 1:nr_sensors
    [~, ~, sensor_noise_result] = sensor_noise(L);
    B.sensor_noise(i, :) = sensor_noise_result;
end

figure;
[pxx_signal, f] = pwelch(B.sensor_noise', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of Sensor noise');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;

%% P2P
biggestP2P = 0;
index = 0;
for i = 1:size(B.brain_signal, 1)
    P2P = max(B.brain_signal(i,:)) - min(B.brain_signal(i,:));
    if P2P > biggestP2P
        biggestP2P = P2P;
        index = i;
    end
end

%% External headmodel
cfg = [];
cfg.method = 'infinite';

headmodel_out = ft_prepare_headmodel(cfg);

%% External sourcemodel
ecg_pos = [0 0 -50];
env_1m_pos = [0 86.6 -50];
env_2m_pos = [-200 0 0];
env_3m_pos = [200 200 -100];
env_10m_pos = [0 1000 0];
pos_out = [ecg_pos; env_1m_pos; env_2m_pos; env_3m_pos; env_10m_pos];

cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = pos_out;

sourcemodel_out = ft_prepare_sourcemodel(cfg);

%% External leadfield
cfg = [];
cfg.grid           = sourcemodel_out;
cfg.headmodel      = headmodel_out;
cfg.grad           = grad;

leadfield_out_prepared = ft_prepare_leadfield(cfg); % compute leadfield
leadfield_out_prepared.leadfield = leadfield_out_prepared.leadfield';
n_sources = length(leadfield_out_prepared.leadfield);
n_sensors = length(leadfield_out_prepared.label);

% Current directions
ecg_vnn =       [1 0 0]; % Outwards
% Environment magnetic field direction perpendicular to vector direction
env_1m_vnn =    [0 1 0]; 
env_2m_vnn =    [1 0 0];
env_3m_vnn =    [1 1 0]/norm([1 1 0]);
env_10m_vnn =   [0 1 0];
vnn_out = [ecg_vnn; env_1m_vnn; env_2m_vnn; env_3m_vnn; env_10m_vnn];
%%
lf_out = zeros(n_sensors, n_sources);
for i = 1:n_sources
    lf_out(:,i) = leadfield_out_prepared.leadfield{i} * vnn_out(i,:)';
end

%% Biological noise
addpath('ECGSimulation')
Q_ecg = 10e-14;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
biological_noise = ecgSimulation(t);

q = zeros(n_sources, n_samples);
q(1, :) = biological_noise;
B.biological_noise = lf_out * q;

figure;

subplot(3,1,1);
plot(t,biological_noise);
title('ECG');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, B.biological_noise)
title('Magnetic field of Biological Noise');
xlabel('Time (s)');
ylabel('T (Tesla)');
grid on;

t_trial = -0.5:(1/Fs):(99.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
biological_noise = ecgSimulation(t);

q = zeros(n_sources, 100000);
q(1, :) = biological_noise;
B.biological_noise = lf_out * q;

subplot(3,1,3);
[pxx_signal, f] = pwelch(B.biological_noise', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of Biological Noise');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;


%% Environmental noise
% Aim for 10e-12 Tesla
% 10-40 Hz

%% 1m
Q_noise = 1e1;
f = 10;
t_trial = -0.5:(1/Fs):(0.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

environmental_noise_result = environmental_noise(t, f, Q_noise);

q = zeros(n_sources, 1000);
q(2, :) = environmental_noise_result;
B.env_1m = lf_out * q;

figure;
subplot(3,1,1);
plot(t, environmental_noise_result);
xlabel('Time (s)');
ylabel('Amplitude');
title('Environmental Noise 1m');
grid on;

subplot(3,1,2);
plot(t, B.env_1m)
title('Magnetic field of environmental noise 1m');
xlabel('Time (s)');
ylabel('T (Tesla)');
grid on;

t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

environmental_noise_result = environmental_noise(t, f, Q_noise);

q = zeros(n_sources, n_samples);
q(2, :) = environmental_noise_result;
B.env_1m = lf_out * q;

subplot(3,1,3);
[pxx_signal, f] = pwelch(B.env_1m', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of environmental noise 1m');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;

%% 2m
Q_noise = 1e2;
f = 16;
t_trial = -0.5:(1/Fs):(0.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
environmental_noise_result = environmental_noise(t, f, Q_noise);

figure;
subplot(3,1,1);
plot(t, environmental_noise_result);
xlabel('Time (s)');
ylabel('Amplitude');
title('Environmental Noise 2m');
grid on;

q = zeros(n_sources, 1000);
q(3, :) = environmental_noise_result;
B.env_2m = lf_out * q;

subplot(3,1,2);
plot(t, B.env_2m)
title('Magnetic field of environmental noise 2m');
xlabel('Time (s)');
ylabel('T (Tesla)');
grid on;

t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
environmental_noise_result = environmental_noise(t, f, Q_noise);

q = zeros(n_sources, n_samples);
q(3, :) = environmental_noise_result;
B.env_2m = lf_out * q;

subplot(3,1,3);
[pxx_signal, f] = pwelch(B.env_2m', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of environmental noise 2m');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;

%% 3m
Q_noise = 1e2;
f = 28;
t_trial = -0.5:(1/Fs):(0.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
environmental_noise_result = environmental_noise(t, f, Q_noise);

figure;
subplot(3,1,1);
plot(t, environmental_noise_result);
xlabel('Time (s)');
ylabel('Amplitude');
title('Environmental Noise 3m');
grid on;

q = zeros(n_sources, 1000);
q(4, :) = environmental_noise_result;
B.env_3m = lf_out * q;

subplot(3,1,2);
plot(t, B.env_3m)
title('Magnetic field of environmental noise 3m');
xlabel('Time (s)');
ylabel('T (Tesla)');
grid on;

t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
environmental_noise_result = environmental_noise(t, f, Q_noise);

q = zeros(n_sources, n_samples);
q(4, :) = environmental_noise_result;
B.env_3m = lf_out * q;

subplot(3,1,3);
[pxx_signal, f] = pwelch(B.env_3m', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of environmental noise 3m');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;

%% 10m
Q_noise = 1e4;
f = 40;
t_trial = -0.5:(1/Fs):(0.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
environmental_noise_result = environmental_noise(t, f, Q_noise);

figure;
subplot(3,1,1);
plot(t, environmental_noise_result);
xlabel('Time (s)');
ylabel('Amplitude');
title('Environmental Noise 10m');
grid on;

q = zeros(n_sources, 1000);
q(5, :) = environmental_noise_result;
B.env_10m = lf_out * q;

subplot(3,1,2);
plot(t, B.env_10m)
title('Magnetic field of environmental noise 10m');
xlabel('Time (s)');
ylabel('T (Tesla)');
grid on;

t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
environmental_noise_result = environmental_noise(t, f, Q_noise);

q = zeros(n_sources, n_samples);
q(5, :) = environmental_noise_result;
B.env_10m = lf_out * q;

subplot(3,1,3);
[pxx_signal, f] = pwelch(B.env_10m', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of environmental noise 10m');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;

%% Total signal
signal = B.brain_signal + B.sensor_noise + B.biological_noise + B.env_1m + B.env_2m + B.env_3m + B.env_10m;

figure;
[pxx_signal, f] = pwelch(signal', [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of signal');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;