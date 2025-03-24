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
data_sss = cell(length(signals), 1); % signals x 1
for i = 1:length(signals)
    cfg = [];
    data_sss{i} = ft_denoise_sss(cfg, data_cell{i});
end

%% Mean Result of all signals
data_sss_mean = data_sss{1}.trial{1};
% data_sss_mean = zeros(length(lf_brain(:,1)), 10000);
% for i = 1:length(signals)
%     data_sss_mean = data_sss_mean + data_sss{i}.trial{1};
% end
% data_sss_mean = data_sss_mean / length(signals);

%% Mean of all signals

data_mean = zeros(length(lf_brain(:,1)), 10000);
for i = 1:length(signals)
    data_mean = data_mean + data_cell{i}.trial{1};
end
data_mean = data_mean / length(signals);

%% Before and after SSS plots
figure;
subplot(2, 1, 1)
[pxx_signal, f] = pwelch(data_mean', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal before SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 1, 2)
[pxx_signal, f] = pwelch(data_sss_mean', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% PSD difference for each component (largest value)
[max_diff, max_one_channel] = psd_diff(data_mean, {data_sss_mean}, freq_dict, Fs, 1);

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
for i = 1:1
    h(i).CData = max_diff(:, i)';
end
xlabel('Sources')
ylabel('Shielding factor (dB)')
grid on
colorbar

figure
h = bar(keys(freq_dict), max_one_channel);
set(h, 'FaceColor', 'flat');
for i = 1:1
    h(i).CData = max_one_channel(:, i)';
end
xlabel('Sources')
ylabel('Shielding factor (dB)')
grid on
colorbar