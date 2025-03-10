% Parameters
Fs = 1000;              % Sampling frequency (Hz)
L = 10000;              % Length of signal (number of samples)
T = 1/Fs;               % Sampling period (s)
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
n_samples = length(t);


%% Sensor noise
[white_noise, pink_noise, sensor_noise_result] = sensor_noise(L);

figure;

subplot(2,1,1);
[pxx_white, f] = pwelch(white_noise, [], [], [], Fs); % Compute PSD
[pxx_pink, ~] = pwelch(pink_noise, [], [], [], Fs); % Compute PSD
pxx_white_T = sqrt(pxx_white);  % Convert to T/Hz^(1/2)
pxx_pink_T = sqrt(pxx_pink);  % Convert to T/Hz^(1/2)
loglog(f, pxx_white_T, 'b');
hold on
loglog(f, pxx_pink_T, 'r');
title('Power Spectral Density of White and Pink Noise');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
legend('White Noise', 'Pink Noise');
grid on;
hold off;

subplot(2,1,2);
[pxx, f] = pwelch(sensor_noise_result, [], [], [], Fs); % Compute PSD
pxx_T = sqrt(pxx);  % Convert to T/Hz^(1/2)
loglog(f, pxx_T);
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
title('Power Spectral Density of Sensor Noise');
grid on;

%% Brain signal

source = 9304; % 1304
t_trial = -0.5:(1/Fs):(0.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
q_signal = brain_signal(source, t);

figure;

subplot(2,1,1);
plot(t, q_signal)
title('Brain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+0.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
q_signal = brain_signal(source, t);

subplot(2,1,2);
[pxx_brain, f] = pwelch(q_signal, [], [], [], Fs); % Compute PSD
pxx_brain_T = sqrt(pxx_brain);  % Convert to T/Hz^(1/2)
loglog(f, pxx_brain_T, 'b');
title('Power Spectral Density of Brain Signal');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;

%% Environmental noise - far
environmental_noise = environmental_noise(t, 10);

figure;
plot(t, environmental_noise);
xlabel('Time (s)');
ylabel('Amplitude');
title('Environmental Noise');
grid on;

%% Biological noise (ecg)
Q_ecg = 10e-14;
biological_noise = Q_ecg*ecgSimulation(t);

figure;

subplot(2,1,1);
plot(t,biological_noise);
title('ECG');
xlabel('Time');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
[pxx_ecg, f] = pwelch(biological_noise, [], [], [], Fs); % Compute PSD
pxx_ecg_T = sqrt(pxx_ecg);  % Convert to T/Hz^(1/2)
loglog(f, pxx_ecg_T, 'b');
title('Power Spectral Density of ECG');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;

%% Signal
signal = sensor_noise_result + q_signal + environmental_noise + biological_noise;
% signal1 = sensor_noise + brain_signal;
% signal2 = signal1 + environmental_noise;
% signal3 = signal2 + biological_noise;

figure;
[pxx_signal, f] = pwelch(signal, [], [], [], Fs); % Compute PSD
pxx_signal_T = sqrt(pxx_signal);  % Convert to T/Hz^(1/2)
loglog(f, pxx_signal_T, 'b');
title('Power Spectral Density of Signal');
xlabel('Frequency (Hz)');
ylabel('PSD (T/Hz^{1/2})');
grid on;