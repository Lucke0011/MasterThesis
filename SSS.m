%% Generate signals
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

[signals, lf_brain, ~, freq_dict] = generate_signals();

%% Data
data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:length(signals)
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
end
%% SSS
data_sss = cell(length(signals), 1); % signals x 1 (5x1)
for i = 1:length(signals)
    cfg = [];
    data_sss{i} = ft_denoise_sss(cfg, data_cell{i});
end

%% Mean of all signals

data_sss_mean = zeros(length(lf_brain(:,1)), 10000);
for i = 1:length(signals)
    data_sss_mean = data_sss_mean + data_sss{i}.trial{1};
end
data_sss_mean = data_sss_mean / length(signals);

%% Before and after SSS plots
figure;
subplot(2, 1, 1)
[pxx_signal, f] = pwelch(data.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal before SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 1, 2)
[pxx_signal, f] = pwelch(data_sss_mean', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal, 'b');
title('PSD of Signal after SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% PSD difference for each component (largest value)
[max_diff, max_one_channel] = psd_diff(data.trial{1}, {data_sss_mean}, freq_dict, Fs, 1);

%% Bar plot of shielding factors

figure
bar(keys(freq_dict), max_diff)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - after of SSS')
grid on

figure
bar(keys(freq_dict), max_one_channel)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - same channel after of SSS')
grid on