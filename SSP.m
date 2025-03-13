%% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

%% Generate signals

[signal, ~, ~, freq_dict] = generate_signals(1);
signal_er = empty_room_signals(1);

% Data
[data, timelock] = generate_data(signal, t);
[data_er, timelock_er] = generate_data(signal_er, t);

signal_fif = 'signal-raw.fif';
signal_er_fif  = 'signal-er-raw.fif';
fieldtrip2fiff(signal_fif, data)
fieldtrip2fiff(signal_er_fif, data_er)

%% Import SSP signal from mne python

fiff_file = 'mne-signal-raw.fif';

cfg         = [];
cfg.dataset = fiff_file;
data_mp     = ft_preprocessing(cfg);
ft_datatype(data_mp)  % should return 'raw'

data_ssp = generate_data(data_mp.trial, t);

%% Before and after SSP plots
figure;
[pxx_signal, f] = pwelch(data.trial{1}', [], [], 0:0.2:500, Fs); % Compute PSD
%pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2) Amplitude Spectral density
loglog(f, pxx_signal, 'b');
title('PSD of Signal before SSP (1 projector)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

figure;
[pxx_signal, f] = pwelch(data_ssp.trial{1}', [], [], 0:0.2:500, Fs); % Compute PSD
%pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2) Amplitude Spectral density
loglog(f, pxx_signal, 'b');
title('PSD of Signal after SSP (1 projector)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% PSD difference for each component (largest value)
[max_diff, max_one_channel] = psd_diff(data.trial{1}, data_ssp.trial{1}, freq_dict, Fs);

%% Box plot of shielding factors

y = [];
for i = 1:length(keys(freq_dict))
    y(i) = max_diff{i}(1,5);
end

figure
bar(keys(freq_dict), y)
xlabel('Source')
ylabel('Shielding factor')
title('Max Shielding factor before/after of SSP (1 projector)')
grid on

y = [];
for i = 1:length(keys(freq_dict))
    y(i) = max_one_channel{i}(1,5);
end

figure
bar(keys(freq_dict), y)
xlabel('Source')
ylabel('Shielding factor')
title('Max Shielding factor before / same channel after of SSP (1 projector)')
grid on

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
