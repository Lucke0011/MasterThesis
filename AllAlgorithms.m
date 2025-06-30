%% Init

% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

%% Generate data
[signals, lf_brain, ~, freq_dict] = generate_signals();
signal_er = empty_room_signals(1);
n_sensors = length(lf_brain(:,1));
n_signals = length(signals);

signal_fif = 'signal-raw-';
signal_er_fif  = 'signal-er-raw';
file_ext = '.fif';

data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:n_signals
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
    signal_file_name = sprintf('%s%d%s', signal_fif, i, file_ext);
    fieldtrip2fiff(signal_file_name, data) % Create .fif files for MNE-Python SSP
end

[data_er, timelock_er] = generate_data(signal_er{1}, t);
empty_room_file_name = sprintf('%s%s', signal_er_fif, file_ext);
fieldtrip2fiff(empty_room_file_name, data_er)

%% HFC
orders = 3;
data_hfc = cell(length(signals), orders); % signals x order
for j = 1:n_signals
    for i = 1:orders
        cfg = [];
        cfg.order = i;
        data_hfc{j, i} = ft_denoise_hfc(cfg, data_cell{j});
    end
end

%% Import SSP signal from mne python
projectors = 5;
folder_path = '/Users/lucke/Exjobb/MNEPython/';
file_prefix = 'mne-signal-raw-';
file_ext = '.fif';
dash = '-';

data_ssp = cell(length(signals), projectors);
for i = 1:n_signals
    for j = 1:projectors
        fileName = sprintf('%s%s%d%s%d%s', folder_path, file_prefix, i, dash, j, file_ext);
        cfg = [];
        cfg.dataset = fileName;
        data_ssp{i, j} = ft_preprocessing(cfg);
        fprintf('Loaded %s\n', fileName);
    end
end

%% Iterative SSS

data_sss = cell(length(signals), 1);
for i = 1:n_signals
    cfg = [];
    cfg.sss.order_in = 9;
    cfg.sss.order_out = 3;
    data_sss{i} = ft_denoise_sss(cfg, data_cell{i});
end

%% AMM

data_amm = cell(n_signals, 1);
for i = 1:n_signals
    cfg = [];
    cfg.channel = 'all';
    data_amm{i} = ft_denoise_amm(cfg, data_cell{i});
end

%% PSD plots

% HFC
figure;
subplot(2, 2, 1)
[psd, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal Before HFC');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 2)
[psd, f] = pwelch(data_hfc{1, 1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal After HFC (1st Order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 3)
[psd, f] = pwelch(data_hfc{1, 2}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal After HFC (2nd Order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 4)
[psd, f] = pwelch(data_hfc{1, 3}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, psd);
title('Signal After HFC (3rd Order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

% SSP
figure;
subplot(2, 2, 1)
[pxx_signal, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal Before SSP');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 2)
[pxx_signal, f] = pwelch(data_ssp{1,1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal After SSP (1 Projector)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 3)
[pxx_signal, f] = pwelch(data_ssp{1,2}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal After SSP (2 Projectors)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 4)
[pxx_signal, f] = pwelch(data_ssp{1,4}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal After SSP (4 Projectors)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

% Iterative SSS
figure;
subplot(2, 1, 1)
[pxx_signal, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal Before Iterative SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 1, 2)
[pxx_signal, f] = pwelch(data_sss{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal After Iterative SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

% AMM
figure;
subplot(2, 1, 1)
[pxx_signal, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal Before AMM');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 1, 2)
[pxx_signal, f] = pwelch(data_amm{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal After AMM');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% Shielding factors

% HFC
hfc_max_diff = zeros(length(keys(freq_dict)), orders);
hfc_max_same_channel = zeros(length(keys(freq_dict)), orders);

i_freqs = [101, 41, 7, 13, 19, 61, 81, 141, 201]; % 20, 8, 1.2, 2.4, 3.6, 12, 26, 28, 40 Hz
n_freqs = length(i_freqs);

for signal = 1:n_signals

    % Before
    [psd_before, ~] = pwelch(data_cell{signal}.trial{1}', [], [], 0:0.2:500, Fs); % 0.2 Hz
    for i_freq = 1:n_freqs
        [psd_before_max, i_sensor] = max(psd_before(i_freqs(i_freq),:));

        % After
        for order = 1:orders
            [psd_after, ~] = pwelch(data_hfc{signal, order}.trial{1}', [], [], 0:0.2:500, Fs);
            [psd_after_max, ~] = max(psd_after(i_freqs(i_freq),:));

            % Shielding factor (max / max)
            psd_diff = psd_before_max / psd_after_max;
    
            % Shielding factor (max / same channel)
            psd_diff_same_channel = psd_before_max / psd_after(i_freqs(i_freq), i_sensor); % Same channel (index) as max
            
            hfc_max_diff(i_freq, order)         = hfc_max_diff(i_freq, order) + psd_diff;
            hfc_max_same_channel(i_freq, order) = hfc_max_same_channel(i_freq, order) + psd_diff_same_channel;
        end
    end
end

ecg_components = 3;
for i = 1:ecg_components-1
    hfc_max_diff(3,:) = hfc_max_diff(3,:) + hfc_max_diff(3+1,:);
    hfc_max_same_channel(3,:) = hfc_max_same_channel(3,:) + hfc_max_same_channel(3+1,:);
    hfc_max_diff(3+1,:) = [];
    hfc_max_same_channel(3+1,:) = [];
end
hfc_max_diff(3,:) = hfc_max_diff(3,:) / ecg_components;
hfc_max_same_channel(3,:) = hfc_max_same_channel(3,:) / ecg_components;

% SSP
ssp_max_diff = zeros(length(keys(freq_dict)), projectors);
ssp_max_same_channel = zeros(length(keys(freq_dict)), projectors);

i_freqs = [101, 41, 7, 13, 19, 61, 81, 141, 201]; % 20, 8, 1.2, 2.4, 3.6, 12, 26, 28, 40 Hz
n_freqs = length(i_freqs);

for signal = 1:n_signals

    % Before
    [psd_before, ~] = pwelch(data_cell{signal}.trial{1}', [], [], 0:0.2:500, Fs); % 0.2 Hz
    for i_freq = 1:n_freqs
        [psd_before_max, i_sensor] = max(psd_before(i_freqs(i_freq),:));

        % After
        for projector = 1:projectors
            [psd_after, ~] = pwelch(data_ssp{signal, projector}.trial{1}', [], [], 0:0.2:500, Fs);
            [psd_after_max, ~] = max(psd_after(i_freqs(i_freq),:));

            % Shielding factor (max / max)
            psd_diff = psd_before_max / psd_after_max;
    
            % Shielding factor (max / same channel)
            psd_diff_same_channel = psd_before_max / psd_after(i_freqs(i_freq), i_sensor); % Same channel (index) as max
            
            ssp_max_diff(i_freq, projector)         = ssp_max_diff(i_freq, projector) + psd_diff;
            ssp_max_same_channel(i_freq, projector) = ssp_max_same_channel(i_freq, projector) + psd_diff_same_channel;
        end
    end
end

ecg_components = 3;
for i = 1:ecg_components-1
    ssp_max_diff(3,:) = ssp_max_diff(3,:) + ssp_max_diff(3+1,:);
    ssp_max_same_channel(3,:) = ssp_max_same_channel(3,:) + ssp_max_same_channel(3+1,:);
    ssp_max_diff(3+1,:) = [];
    ssp_max_same_channel(3+1,:) = [];
end
ssp_max_diff(3,:) = ssp_max_diff(3,:) / ecg_components;
ssp_max_same_channel(3,:) = ssp_max_same_channel(3,:) / ecg_components;

% Iterative SSS
sss_max_diff = zeros(length(keys(freq_dict)), 1);
sss_max_same_channel = zeros(length(keys(freq_dict)), 1);

i_freqs = [101, 41, 7, 13, 19, 61, 81, 141, 201]; % 20, 8, 1.2, 2.4, 3.6, 12, 26, 28, 40 Hz
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
            
            sss_max_diff(i_freq, order)         = sss_max_diff(i_freq, order) + psd_diff;
            sss_max_same_channel(i_freq, order) = sss_max_same_channel(i_freq, order) + psd_diff_same_channel;
        end
    end
end

ecg_components = 3;
for i = 1:ecg_components-1
    sss_max_diff(3,:) = sss_max_diff(3,:) + sss_max_diff(3+1,:);
    sss_max_same_channel(3,:) = sss_max_same_channel(3,:) + sss_max_same_channel(3+1,:);
    sss_max_diff(3+1,:) = [];
    sss_max_same_channel(3+1,:) = [];
end
sss_max_diff(3,:) = sss_max_diff(3,:) / ecg_components;
sss_max_same_channel(3,:) = sss_max_same_channel(3,:) / ecg_components;

% AMM
amm_max_diff = zeros(length(keys(freq_dict)), 1);
amm_max_same_channel = zeros(length(keys(freq_dict)), 1);

i_freqs = [101, 41, 7, 13, 19, 61, 81, 141, 201]; % 20, 8, 1.2, 2.4, 3.6, 12, 26, 28, 40 Hz
n_freqs = length(i_freqs);

for signal = 1:n_signals

    % Before
    [psd_before, ~] = pwelch(data_cell{signal}.trial{1}', [], [], 0:0.2:500, Fs); % 0.2 Hz
    for i_freq = 1:n_freqs
        [psd_before_max, i_sensor] = max(psd_before(i_freqs(i_freq),:));

        % After
        for order = 1:1
            [psd_after, ~] = pwelch(data_amm{signal, order}.trial{1}', [], [], 0:0.2:500, Fs);
            [psd_after_max, ~] = max(psd_after(i_freqs(i_freq),:));

            % Shielding factor (max / max)
            psd_diff = psd_before_max / psd_after_max;
    
            % Shielding factor (max / same channel)
            psd_diff_same_channel = psd_before_max / psd_after(i_freqs(i_freq), i_sensor); % Same channel (index) as max
            
            amm_max_diff(i_freq, order)         = amm_max_diff(i_freq, order) + psd_diff;
            amm_max_same_channel(i_freq, order) = amm_max_same_channel(i_freq, order) + psd_diff_same_channel;
        end
    end
end

ecg_components = 3;
for i = 1:ecg_components-1
    amm_max_diff(3,:) = amm_max_diff(3,:) + amm_max_diff(3+1,:);
    amm_max_same_channel(3,:) = amm_max_same_channel(3,:) + amm_max_same_channel(3+1,:);
    amm_max_diff(3+1,:) = [];
    amm_max_same_channel(3+1,:) = [];
end
amm_max_diff(3,:) = amm_max_diff(3,:) / ecg_components;
amm_max_same_channel(3,:) = amm_max_same_channel(3,:) / ecg_components;

% Remove ecg 2 and 3
freq_dict = remove(freq_dict, "ecg 2");
freq_dict = remove(freq_dict, "ecg 3");

% Mean and to dB
% HFC
hfc_max_diff = hfc_max_diff / n_signals;
hfc_max_diff = 10*log10(hfc_max_diff);
hfc_max_same_channel = hfc_max_same_channel / n_signals;
hfc_max_same_channel = 10*log10(hfc_max_same_channel);

% SSP
ssp_max_diff = ssp_max_diff / n_signals;
ssp_max_diff = 10*log10(ssp_max_diff);
ssp_max_same_channel = ssp_max_same_channel / n_signals;
ssp_max_same_channel = 10*log10(ssp_max_same_channel);

% Iterative SSS
sss_max_diff = sss_max_diff / n_signals;
sss_max_diff = 10*log10(sss_max_diff);
sss_max_same_channel = sss_max_same_channel / n_signals;
sss_max_same_channel = 10*log10(sss_max_same_channel);

% AMM
amm_max_diff = amm_max_diff / n_signals;
amm_max_diff = 10*log10(amm_max_diff);
amm_max_same_channel = amm_max_same_channel / n_signals;
amm_max_same_channel = 10*log10(amm_max_same_channel);

%% Bar plot of shielding factors

% HFC
figure
bar(keys(freq_dict), hfc_max_diff)
title('Shielding Factors of HFC');
xlabel('Source')
ylabel('Shielding factor (dB)')
grid on

% SSP
figure
bar(keys(freq_dict), ssp_max_diff)
title('Shielding Factors of SSP');
xlabel('Source')
ylabel('Shielding factor (dB)')
grid on

% Iterative SSS
figure
bar(keys(freq_dict), sss_max_diff)
title('Shielding Factors of Iterative SSS');
xlabel('Source')
ylabel('Shielding factor (dB)')
grid on

% AMM
figure
bar(keys(freq_dict), amm_max_diff)
title('Shielding Factors of AMM');
xlabel('Source')
ylabel('Shielding factor (dB)')
grid on

%% Brain signal shielding factors per algorithm

algorithms = {'1st order HFC', '2nd order HFC', '3rd order HFC', '4 Projectors SSP', 'Iterative SSS', 'AMM'};
shieldingFactors = [hfc_max_diff(1, 1), hfc_max_diff(1, 2), hfc_max_diff(1, 3), ssp_max_diff(1, 4), sss_max_diff(1), amm_max_diff(1)];

figure
bar(algorithms, shieldingFactors)
title('Brain Signal Shielding Factors per Algorithm');
xlabel('Algorithms')
ylabel('Shielding factor (dB)')
grid on