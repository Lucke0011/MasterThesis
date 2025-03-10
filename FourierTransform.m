%% Fourier transform to get components
% Parameters
Fs = 1000;          % Sampling frequency (Hz)
L = 10000;          % Length of signal
components = 3;

[ecg_ft_result, locs] = ecg_ft(components, Fs, L);

figure
[pxx_signal, f] = pwelch(ecg_ft_result, [], [], [], Fs); % Compute PSD
% pxx_signal_T = sqrt(pxx_signal); % Convert to T/Hz^(1/2) Amplitude Spectral density
loglog(f, pxx_signal_T, 'b');
title('PSD of reconstructed ecg');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
