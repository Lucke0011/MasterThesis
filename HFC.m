%% Init
load('grad.mat')

% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

[signals, lf_brain, ~, freq_dict] = generate_signals();
n_sensors = length(lf_brain(:,1));
n_signals = length(signals);

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
for i = 1:n_signals
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
end

%% Mean of raw signals

data_mean = zeros(n_sensors, 10000);
for i = 1:n_signals
    data_mean = data_mean + data_cell{i}.trial{1};
end
data_mean = data_mean / length(signals);

%% HFC
orders = 3;
data_hfc_cell = cell(length(signals), orders); % signals x order
for j = 1:n_signals
    for i = 1:orders
        cfg = [];
        cfg.order = i;
        data_hfc_cell{j, i} = ft_denoise_hfc(cfg, data_cell{j});
    end
end

%% Before and after HFC plots
figure;
subplot(2, 2, 1)
[psd, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal before HFC');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 2)
[psd, f] = pwelch(data_hfc_cell{1, 1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal after HFC (1st order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 3)
[psd, f] = pwelch(data_hfc_cell{1, 2}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal after HFC (2nd order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 4)
[psd, f] = pwelch(data_hfc_cell{1, 3}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal after HFC (3rd order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% Mean of Shielding factors 1 min 9 sec (old was 1 min 29 sec)

max_diff = zeros(length(keys(freq_dict)), orders);
max_same_channel = zeros(length(keys(freq_dict)), orders);

n_freqs = length(i_freqs);
i_freqs = [61, 81, 141, 201, 7, 13, 19, 101]; % 12, 26, 28, 40, 1.2, 2.4, 3.6, 20 Hz

for signal = 1:n_signals

    % Before
    [psd_before, ~] = pwelch(data_cell{signal}.trial{1}', [], [], 0:0.2:500, Fs); % 0.2 Hz
    for i_freq = 1:n_freqs
        [psd_before_max, i_sensor] = max(psd_before(i_freqs(i_freq),:));

        % After
        for order = 1:orders
            [psd_after, ~] = pwelch(data_hfc_cell{signal, order}.trial{1}', [], [], 0:0.2:500, Fs);
            [psd_after_max, ~] = max(psd_after(i_freqs(i_freq),:));

            % Shielding factor (max / max)
            psd_diff = psd_before_max / psd_after_max;
    
            % Shielding factor (max / same channel)
            psd_diff_same_channel = psd_before_max / psd_after(i_freqs(i_freq), i_sensor); % Same channel (index) as max
            
            max_diff(i_freq, order)         = max_diff(i_freq, order) + psd_diff;
            max_same_channel(i_freq, order) = max_same_channel(i_freq, order) + psd_diff_same_channel;
        end
    end
end

ecg_components = 3;
for i = 1:ecg_components-1
    max_diff(5,:) = max_diff(5,:) + max_diff(5+1,:);
    max_same_channel(5,:) = max_same_channel(5,:) + max_same_channel(5+1,:);
    max_diff(5+1,:) = [];
    max_same_channel(5+1,:) = [];
end
max_diff(5,:) = max_diff(5,:) / ecg_components;
max_same_channel(5,:) = max_same_channel(5,:) / ecg_components;

% Remove ecg 2 and 3
freq_dict = remove(freq_dict, "ecg 2");
freq_dict = remove(freq_dict, "ecg 3");
freq_dict("ecg") = freq_dict("ecg 1");
freq_dict("brain signal") = freq_dict("brain_signal");
freq_dict = remove(freq_dict, "brain_signal");
freq_dict = remove(freq_dict, "ecg 1");

max_diff = max_diff / length(signals);
max_diff = 20*log10(max_diff);
max_same_channel = max_same_channel / length(signals);
max_same_channel = 20*log10(max_same_channel);

%% Bar plot of shielding factors

figure
bar(keys(freq_dict), max_diff)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - after of HFC')
grid on

figure
bar(keys(freq_dict), max_same_channel)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - same channel after of HFC')
grid on

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