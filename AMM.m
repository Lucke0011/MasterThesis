%% Generate signals
Fs = 1000;
t_trial = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
t = (0:(length(t_trial)-1)) * (1/Fs);

[signals, lf_brain, ~, freq_dict] = generate_signals(); % Hitting noise floor on 10m 
n_sensors = length(lf_brain(:,1));

%% Data
data_cell = cell(length(signals), 1);
timelock_cell = cell(length(signals), 1);
for i = 1:length(signals)
    [data, timelock] = generate_data(signals{i}, t);
    data_cell{i} = data;
    timelock_cell{i} = timelock;
end

%% AMM spm

data_amm = cell(length(signals), 1);
for i = 1:length(signals)
    D = spm_eeg_ft2spm(data_cell{i}, 'spm_raw');
    S = [];
    S.D = D;
    [mfD,Yinds] = spm_opm_amm(S);
end

%% AMM fieldtrip

data_amm = cell(length(signals), 1);
for i = 1:length(signals)
    cfg = [];
    cfg.channel = '*bz';
    data_amm{i} = ft_denoise_amm(cfg, data_cell{i});
end

%% Before and after SSS plots
figure;
subplot(2, 1, 1)
[pxx_signal, f] = pwelch(data_cell{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal before AMM');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

subplot(2, 1, 2)
[pxx_signal, f] = pwelch(data_amm{1}.trial{1}', [], [], 0:0.2:500, Fs);
loglog(f, pxx_signal);
title('Signal after AMM');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;

%% Mean of Shielding factors
max_diff = zeros(length(keys(freq_dict)), 1);
max_one_channel = zeros(length(keys(freq_dict)), 1);

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
            temp = interp1(f, psd(:,j), freq);
            if temp > psd_before
                psd_before = temp;
                index = j;
            end
        end
            
        % After
        for order = 1:1
            [psd, f] = pwelch(data_amm{signal, order}.trial{1}', [], [], 0:0.2:500, Fs);
            
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
            
            max_diff(i, order)        = max_diff(i, order) + psd_diff;
            max_one_channel(i, order) = max_one_channel(i, order) + psd_diff_one_channel;
        end
    end
end

% mean of ecg components
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
bar(keys(freq_dict), max_diff)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - after of AMM')
grid on

figure
bar(keys(freq_dict), max_one_channel)
xlabel('Source')
ylabel('Shielding factor (dB)')
title('Max Shielding factor before - same channel after of AMM')
grid on