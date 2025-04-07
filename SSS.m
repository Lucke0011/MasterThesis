%% Generate signals
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

[signals, lf_brain, ~, freq_dict] = generate_signals();
n_sensors = length(lf_brain(:,1));
n_signals = length(signals);

%% Data
data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:length(signals)
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
end

%% MNE-SSS
signal_fif = 'sss-signal-raw-';
file_ext = '.fif';

data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:n_signals
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
    signal_file_name = sprintf('%s%d%s', signal_fif, i, file_ext);
    fieldtrip2fiff(signal_file_name, data, 'hdr', hdr)
end

%% Import SSS from MNE

folder_path = '/Users/lucke/Exjobb/MNEPython/';
file_prefix = 'mne-sss-signal-raw-';
file_ext = '.fif';

data_sss = cell(length(signals), 1);
for i = 1:n_signals
    fileName = sprintf('%s%s%d%s', folder_path, file_prefix, i, file_ext);
    cfg = [];
    cfg.dataset = fileName;
    data_sss{i} = ft_preprocessing(cfg);
end

%% Mean of raw signals

data_mean = zeros(n_sensors, 10000);
for i = 1:n_signals
    data_mean = data_mean + data_cell{i}.trial{1};
end
data_mean = data_mean / length(signals);

%% Fieldtrip SSS

data_sss = cell(length(signals), 1); % signals x 1
for i = 1:n_signals
    cfg = [];
    cfg.sss.order_in = 9;
    cfg.sss.order_out = 3;
    cfg.sss.regularize = 1;
    data_sss{i} = ft_denoise_sss(cfg, data_cell{i});
end

%% Before and after SSS plots
figure;
subplot(2, 1, 1)
[pxx_signal, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal before SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 1, 2)
[pxx_signal, f] = pwelch(data_sss{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% Mean of Shielding factors

max_diff = zeros(length(keys(freq_dict)), 1);
max_same_channel = zeros(length(keys(freq_dict)), 1);

i_freqs = [61, 81, 141, 201, 7, 13, 19, 101]; % 12, 26, 28, 40, 1.2, 2.4, 3.6, 20 Hz
n_freqs = length(i_freqs);

for signal = 1:n_signals

    % Before
    [psd_before, ~] = pwelch(data_cell{signal}.trial{1}', [], [], 0:0.2:500, Fs); % 0.2 Hz
    for i_freq = 1:n_freqs
        [psd_before_max, i_sensor] = max(psd_before(i_freqs(i_freq),:));

        % After
        for order = 1:1
            [psd_after, ~] = pwelch(data_sss{signal, order}.trial{1}', [], [], 0:0.2:500, Fs);
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

max_diff = max_diff / n_signals;
max_diff = 20*log10(max_diff);
max_same_channel = max_same_channel / n_signals;
max_same_channel = 20*log10(max_same_channel);

%% Bar plot of shielding factors

figure
bar(keys(freq_dict), max_diff)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - after of SSS')
grid on

figure
bar(keys(freq_dict), max_same_channel)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - same channel after of SSS')
grid on