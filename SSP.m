%% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
projectors = 10;

%% Generate signals

[signals, lf_brain, ~, freq_dict] = generate_signals(); % Run with normal values for SSP
signal_er = empty_room_signals(1);

%% Data
signal_fif = 'signal-raw-';
signal_er_fif  = 'signal-er-raw';
file_ext = '.fif';

data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:length(signals)
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
    signal_file_name = sprintf('%s%d%s', signal_fif, i, file_ext);
    fieldtrip2fiff(signal_file_name, data)
end

[data_er, timelock_er] = generate_data(signal_er{1}, t);
empty_room_file_name = sprintf('%s%s', signal_er_fif, file_ext);
fieldtrip2fiff(empty_room_file_name, data_er)

%% Import SSP signal from mne python

folder_path = '/Users/lucke/Exjobb/MNEPython/';
file_prefix = 'mne-signal-raw-';
file_ext = '.fif';
dash = '-';

data_ssp = cell(length(signals), projectors);
for i = 1:length(signals)
    for j = 1:projectors
        fileName = sprintf('%s%s%d%s%d%s', folder_path, file_prefix, i, dash, j, file_ext);
        cfg = [];
        cfg.dataset = fileName;
        data_ssp{i, j} = ft_preprocessing(cfg);
        fprintf('Loaded %s\n', fileName);
    end
end

%% Mean of all signals

data_ssp_mean = cell(projectors, 1);
for i = 1:projectors
    mean_temp = zeros(length(lf_brain(:,1)), 10000);
    for j = 1:length(signals)
        mean_temp = mean_temp + data_ssp{j, i}.trial{1};
    end
    data_ssp_mean{i} = mean_temp / length(signals);
end

%% Before and after SSP plots
figure;
subplot(2, 2, 1)
[pxx_signal, f] = pwelch(data.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal before SSP');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 2)
[pxx_signal, f] = pwelch(data_ssp_mean{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal after SSP (1 projector)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 3)
[pxx_signal, f] = pwelch(data_ssp_mean{3}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal after SSP (3 projectors)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 4)
[pxx_signal, f] = pwelch(data_ssp_mean{10}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal after SSP (10 projectors)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% PSD difference for each component (largest value)
[max_diff, max_one_channel] = psd_diff(data.trial{1}, data_ssp_mean, freq_dict, Fs, projectors);

%% Box plot of shielding factors

% y = [];
% for i = 1:length(keys(freq_dict))
%     y(i) = max_diff{i}(1,5);
% end

figure
bar(keys(freq_dict), max_diff)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - after of SSP')
grid on

% y = [];
% for i = 1:length(keys(freq_dict))
%     y(i) = max_one_channel{i}(1,5);
% end

figure
bar(keys(freq_dict), max_one_channel)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - same channel after of SSP')
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
