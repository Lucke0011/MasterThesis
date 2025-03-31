%% Init
% Parameters
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);
projectors = 10;

% Generate signals
[signals, lf_brain, ~, freq_dict] = generate_signals(); % Run with normal values for SSP
n_sensors = length(lf_brain(:,1));
signal_er = empty_room_signals(1);

%% Data
signal_fif = 'signal-raw-';
signal_er_fif  = 'signal-er-raw';
file_ext = '.fif';

data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:length(signals)
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
    signal_file_name = sprintf('%s%d%s', signal_fif, i, file_ext);
    %fieldtrip2fiff(signal_file_name, data)
end

%% t
[data_er, timelock_er] = generate_data(signal_er{1}, t);
empty_room_file_name = sprintf('%s%s', signal_er_fif, file_ext);
fieldtrip2fiff(empty_room_file_name, data_er)

%% Mean of raw signals

data_mean = zeros(n_sensors, 10000);
for i = 1:length(signals)
    data_mean = data_mean + data_cell{i}.trial{1};
end
data_mean = data_mean / length(signals);

%% Import SSP signal from mne python

folder_path = '/Users/lucke/Exjobb/MNEPython/';
file_prefix = 'mne-signal-raw-';
file_ext = '.fif';
dash = '-';

data_ssp = cell(length(signals), projectors);
for i = 1:length(signals)
    for j = 1:projectors
        fileName = sprintf('%s%s%d%s%d%s', folder_path, file_prefix, i, dash, j, file_ext);
        cfg = [];
        cfg.dataset = fileName;
        data_ssp{i, j} = ft_preprocessing(cfg);
        fprintf('Loaded %s\n', fileName);
    end
end

%% Before and after SSP plots
figure;
subplot(2, 2, 1)
[pxx_signal, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal before SSP');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 2)
[pxx_signal, f] = pwelch(data_ssp{1,1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after SSP (1 projector)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 3)
[pxx_signal, f] = pwelch(data_ssp{1,3}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after SSP (3 projectors)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 2, 4)
[pxx_signal, f] = pwelch(data_ssp{1,10}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after SSP (10 projectors)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% Mean of Shielding factors
max_diff = zeros(length(keys(freq_dict)), projectors);
max_one_channel = zeros(length(keys(freq_dict)), projectors);

freqs = keys(freq_dict);
for i = 1:length(freqs)
    freq = freq_dict(freqs(i));

    for signal = 1:length(signals)
        % Before
        [psd, f] = pwelch(data_cell{signal}.trial{1}', [], [], 0:0.2:500, Fs); % 0.2 Hz
        
        % Highest value measured
        psd_before = 0;
        index = 0;
        for j = 1:123
            temp = interp1(f, psd(:,j), freq); % change to max(psd(index of freq,:))
            if temp > psd_before
                psd_before = temp;
                index = j;
            end
        end
            
        % After
        for projector = 1:projectors
            [psd, f] = pwelch(data_ssp{signal, projector}.trial{1}', [], [], 0:0.2:500, Fs);
            
            % Highest value measured
            psd_after = 0;
            for j = 1:n_sensors
                temp = interp1(f, psd(:,j), freq);
                if temp > psd_after
                    psd_after = temp;
                end
            end
            
            % Shielding factor (max / max)
            psd_diff = psd_before / psd_after;
    
            % Shielding factor (max / same channel)
            psd_diff_one_channel = psd_before / interp1(f, psd(:,index), freq); % Same channel (index) as max
            
            max_diff(i, projector)        = max_diff(i, projector) + psd_diff;
            max_one_channel(i, projector) = max_one_channel(i, projector) + psd_diff_one_channel;
        end
    end
end

% Move mean of ecg here
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

max_diff = max_diff / length(signals);
max_diff = 20*log10(max_diff);
max_one_channel = max_one_channel / length(signals);
max_one_channel = 20*log10(max_one_channel);

%% Bar plot of shielding factors

figure
h = bar(keys(freq_dict), max_diff);
set(h, 'FaceColor', 'flat');
for i = 1:projectors
    h(i).CData = max_diff(:, i)';
end
xlabel('Sources')
ylabel('Shielding factor (dB)')
grid on
colormap hot

figure
h = bar(keys(freq_dict), max_one_channel);
set(h, 'FaceColor', 'flat');
for i = 1:projectors
    h(i).CData = max_one_channel(:, i)';
end
xlabel('Sources')
ylabel('Shielding factor (dB)')
grid on
colormap(flipud(hot))

%% g

figure
bar(keys(freq_dict), max_diff)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - after of HFC')
grid on

figure
bar(keys(freq_dict), max_one_channel)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - same channel after of HFC')
grid on

%% Topoplot before and after

load('grad.mat')

%before
cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';

ft_topoplotER(cfg, timelock);

%after
timelock_ssp        = [];
timelock_ssp.dimord = 'chan_time';
timelock_ssp.avg    = data_ssp{1,10}.trial{1};
timelock_ssp.label  = grad.label;
timelock_ssp.time   = t;
timelock_ssp.grad   = grad;

timelock_ssp = ft_datatype_timelock(timelock_ssp);

ft_topoplotER(cfg, timelock_ssp);
