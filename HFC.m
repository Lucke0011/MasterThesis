% clear all

%% Load files
load('grad.mat')

% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
n_trials = 1;
t = (0:(length(t_trial)-1)) * (1/Fs);w

[signals, ~, ~, freq_dict] = generate_signals(n_trials);

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

%% HFC
orders = 3;
data_hfc = cell(orders, 1);
for i = 1:orders
    cfg = [];
    cfg.order = i;
    data_hfc{i} = ft_denoise_hfc(cfg,data);
end

%% Before and after HFC plots
figure;
subplot(2, 2, 1)
[pxx_signal, f] = pwelch(data.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal before HFC');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 2)
[pxx_signal, f] = pwelch(data_hfc{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal after HFC (1st order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 3)
[pxx_signal, f] = pwelch(data_hfc{2}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal after HFC (2nd order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 4)
[pxx_signal, f] = pwelch(data_hfc{3}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal after HFC (3rd order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% PSD difference for each component (largest value)
[max_diff, max_one_channel] = psd_diff(data.trial{1}, data_hfc, freq_dict, Fs, orders);

%% Shielding factor
% hfc = struct();
% 
% hfc_result = cell(length(keys(freq_dict)), 1);
% shielding_factor_max_dict = dictionary();
% shielding_factor_median_dict = dictionary();
% 
% [pxx_signal, f] = pwelch(data.trial{1}', [], [], 0:0.2:500, Fs);
% psd_signal_before = pxx_signal';
% 
% [pxx_signal, f] = pwelch(data_hfc.trial{1}', [], [], 0:0.2:500, Fs);
% psd_signal_after = pxx_signal';
% 
% freqs = keys(freq_dict);
% for i = 1:length(freqs)
%     freq = freq_dict(freqs(i));
% 
%     idx = find(abs(f - freq) < 1e-6); % Adjust tolerance as needed
% 
%     signal_before = psd_signal_before(:,idx);
%     signal_after = psd_signal_after(:,idx);
% 
%     shielding_factor = signal_before ./ signal_after;
% 
%     shielding_factor_max = max(shielding_factor);
%     shielding_factor_median = median(shielding_factor);
% 
%     hfc_result{i} = [signal_before, signal_after, shielding_factor];
%     shielding_factor_max_dict(freqs(i)) = shielding_factor_max;
%     shielding_factor_median_dict(freqs(i)) = shielding_factor_median;
% end
% 
% hfc.hfc_result = hfc_result;
% hfc.shielding_factor_max_dict = shielding_factor_max_dict;
% hfc.shielding_factor_median_dict = shielding_factor_median_dict;

%% Box plot of shielding factors

% y = [];
% for i = 1:length(keys(freq_dict))
%     y(i) = max_diff{i}(1,5);
% end

figure
bar(keys(freq_dict), max_diff)
xlabel('Source')
ylabel('Shielding factor')
title('Max Shielding factor before/after of HFC (3rd order)')
grid on

% y = [];
% for i = 1:length(keys(freq_dict))
%     y(i) = max_one_channel{i}(1,5);
% end

figure
bar(keys(freq_dict), max_one_channel)
xlabel('Source')
ylabel('Shielding factor')
title('Max Shielding factor before / same channel after of HFC (3rd order)')
grid on

% figure
% bar(keys(shielding_factor_median_dict), values(shielding_factor_median_dict))
% xlabel('Source')
% ylabel('Shielding factor')
% title('Median Shielding factor of HFC')
% grid on
% 
% figure
% bar(keys(shielding_factor_max_dict), values(shielding_factor_max_dict))
% xlabel('Source')
% ylabel('Shielding factor')
% title('Max Shielding factor of HFC')
% grid on

%% Topoplot before and after

%before
cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';

ft_topoplotER(cfg, timelock);

%after
timelock_hfc        = [];
timelock_hfc.dimord = 'chan_time';
timelock_hfc.avg    = data_hfc.trial{1};
timelock_hfc.label  = grad.label;
timelock_hfc.time   = t;
timelock_hfc.grad   = grad;

timelock_hfc = ft_datatype_timelock(timelock_hfc);

ft_topoplotER(cfg, timelock_hfc);