Fs = 5000;
Fs_new = 1000;
M = Fs / Fs_new;

%% Preprocessing
cfg = [];
cfg.datafile        = 'EmptyRoomOPM_raw.fif'; % Empty room without participant --> sensor and room noise
cfg.channel         = '*bz';
empty_room_data = ft_preprocessing(cfg);

%%

cfg = [];
cfg.datafile        = 'RSECOPM_raw.fif'; % Resting state, eyes closed --> sensor, room and brain noise
cfg.channel         = '*bz';
resting_state_data = ft_preprocessing(cfg);

%%

n_channels_er = length(empty_room_data.label);
n_samples_er = length(empty_room_data.time{1});
f_sample_er = empty_room_data.fsample;
empty_room_signal = empty_room_data.trial{1}; % n_channels x n_samples
empty_room_signal = decimate(empty_room_signal(76,:), M);

%%
n_channels_rs = length(resting_state_data.label);
n_samples_rs = length(resting_state_data.time{1});
f_sample_rs = resting_state_data.fsample;
resting_state_signal = resting_state_data.trial{1}; % n_channels x n_samples
resting_state_signal = decimate(resting_state_signal(76,:), M);

%% Plots

figure
plot(downsample(empty_room_data.time{1}, M), empty_room_signal)
xlabel('time')
ylabel('Amplitude')
title('Empty Room')
grid on;

%%

figure
plot(resting_state_data.time{1}, resting_state_signal, 'b')
xlabel('time')
ylabel('Amplitude')
title('Resting State')
grid on;

%% PSD Empty Room compared to sensor noise

[~, ~, sensor_noise_result] = sensor_noise(10000);

figure
[pxx, f] = pwelch(empty_room_signal, [], [], [], Fs_new); % Compute PSD
pxx_T = sqrt(pxx);  % Convert to T/Hz^(1/2)
loglog(f, pxx_T, 'b');
hold on;
[pxxn, f] = pwelch(sensor_noise_result, [], [], [], Fs_new); % Compute PSD
pxx_n_T = sqrt(pxxn);  % Convert to T/Hz^(1/2)
loglog(f, pxx_n_T, 'r');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
title('Empty Room compared to Sensor Noise');
legend('Empty Room', 'Sensor Noise');
grid on;

%% PSD Empty Room compared to simulation

figure
[pxx, f] = pwelch(empty_room_signal, [], [], [], Fs_new); % Compute PSD
pxx_T = sqrt(pxx);  % Convert to T/Hz^(1/2)
loglog(f, pxx_T, 'b');
hold on;
[pxxn, f] = pwelch(signal(30,:), [], [], [], Fs_new); % Compute PSD
pxx_n_T = sqrt(pxxn);  % Convert to T/Hz^(1/2)
loglog(f, pxx_n_T, 'r');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
title('Empty Room Signal Compared to Simulated Signal');
legend('Empty Room', 'Simulated');
grid on;