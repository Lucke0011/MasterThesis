% clear all

%% Load files
load('grad.mat')

% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

[signals, lf_brain, ~, freq_dict] = generate_signals();

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
data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:length(signals)
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
end
%% HFC
orders = 3;
data_hfc_cell = cell(length(signals), orders); % signals x order
for j = 1:length(signals)
    for i = 1:orders
        cfg = [];
        cfg.order = i;
        data_hfc_cell{j, i} = ft_denoise_hfc(cfg, data_cell{j});
    end
end

%% Mean result of all signals

data_hfc_mean = cell(orders, 1);
for i = 1:orders
    mean_temp = zeros(length(lf_brain(:,1)), 10000);
    for j = 1:length(signals)
        mean_temp = mean_temp + data_hfc_cell{j, i}.trial{1};
    end
    data_hfc_mean{i} = mean_temp / length(signals);
end

%% Mean of all signals

data_mean = zeros(length(lf_brain(:,1)), 10000);
for i = 1:length(signals)
    data_mean = data_mean + data_cell{i}.trial{1};
end
data_mean = data_mean / length(signals);

%% Before and after HFC plots
figure;
subplot(2, 2, 1)
[pxx_signal, f] = pwelch(data_mean', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal before HFC');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 2)
[pxx_signal, f] = pwelch(data_hfc_mean{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after HFC (1st order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 3)
[pxx_signal, f] = pwelch(data_hfc_mean{2}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after HFC (2nd order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 4)
[pxx_signal, f] = pwelch(data_hfc_mean{3}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after HFC (3rd order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% PSD difference for each component (largest value)
[max_diff, max_one_channel] = psd_diff(data_mean, data_hfc_mean, freq_dict, Fs, orders);

%% Mean of ecg

ecg_components = 3;
for i = 1:ecg_components-1
    max_diff(5,:) = max_diff(5,:) + max_diff(5+1,:);
    max_one_channel(5,:) = max_one_channel(5,:) + max_one_channel(5+1,:);
    max_diff(5+1,:) = [];
    max_one_channel(5+1,:) = [];
end
max_diff(5,:) = max_diff(5,:) / ecg_components;
max_one_channel(5,:) = max_one_channel(5,:) / ecg_components;

% Remove ecg 2 and 3
freq_dict = remove(freq_dict, "ecg 2");
freq_dict = remove(freq_dict, "ecg 3");
freq_dict("ecg") = freq_dict("ecg 1");
freq_dict("brain signal") = freq_dict("brain_signal");
freq_dict = remove(freq_dict, "brain_signal");
freq_dict = remove(freq_dict, "ecg 1");

%% Bar plot of shielding factors

figure
h = bar(keys(freq_dict), max_diff);
set(h, 'FaceColor', 'flat');
for i = 1:orders
    h(i).CData = max_diff(:, i)';
end
xlabel('Sources')
ylabel('Shielding factor (dB)')
grid on
colorbar

figure
h = bar(keys(freq_dict), max_one_channel);
set(h, 'FaceColor', 'flat');
for i = 1:orders
    h(i).CData = max_one_channel(:, i)';
end
xlabel('Sources')
ylabel('Shielding factor (dB)')
grid on
colorbar

%% Topoplot before and after

%before
cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';

ft_topoplotER(cfg, timelock_cell);

%after
timelock_hfc        = [];
timelock_hfc.dimord = 'chan_time';
timelock_hfc.avg    = data_hfc.trial{1};
timelock_hfc.label  = grad.label;
timelock_hfc.time   = t;
timelock_hfc.grad   = grad;

timelock_hfc = ft_datatype_timelock(timelock_hfc);

ft_topoplotER(cfg, timelock_hfc);