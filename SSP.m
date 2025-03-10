% clear all

%% Load files
load('grad.mat')

% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
n_trials = 1;
t = (0:(length(t_trial)-1)) * (1/Fs);

[signals, ~, ~, noise_freq_dict] = generate_signals(n_trials);

%% Butter worth 4th Lowpass filter
% fc = 10;  % Cutoff frequency
% order = 4; % Filter order
% signal_lp = cell(n_trials,1);
% [b, a] = butter(order, fc/(Fs/2), 'high'); % Butterworth filter design
% for i = 1:n_trials
%     for j = 1:123
%         signal_lp{i}(j,:) = filtfilt(b, a, signal{i}(j,:)); % Apply zero-phase filtering
%     end
% end

%% Data
[data, timelock] = generate_data(signals, t);

%% FT Preprocessing
cfg = [];
data_preprocessed = ft_preprocessing(cfg, data);

%% SSP
cfg = [];
cfg.ssp = 'all';
data_ssp = ft_denoise_ssp(cfg,data_preprocessed);

%% Before and after HFC plots
figure;
[pxx_signal, f] = pwelch(data.trial{1}', [], [], 0:0.2:500, Fs); % Compute PSD
%pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2) Amplitude Spectral density
loglog(f, pxx_signal, 'b');
title('Power Spectral Density of Signal before SSP');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

figure;
[pxx_signal, f] = pwelch(data_ssp.trial{1}', [], [], 0:0.2:500, Fs); % Compute PSD
%pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2) Amplitude Spectral density
loglog(f, pxx_signal, 'b');
title('Power Spectral Density of Signal after SSP');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% PSD difference for each component
psd_diff_result = psd_diff(data.trial{1}, data_ssp.trial{1}, noise_freq_dict, Fs);

%% Topoplot before and after

%before
cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';

ft_topoplotER(cfg, timelock);

%after
timelock_ssp        = [];
timelock_ssp.dimord = 'chan_time';
timelock_ssp.avg    = data_ssp.trial{1};
timelock_ssp.label  = grad.label;
timelock_ssp.time   = t;
timelock_ssp.grad   = grad;

timelock_ssp = ft_datatype_timelock(timelock_ssp);

ft_topoplotER(cfg, timelock_ssp);
