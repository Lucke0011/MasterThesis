function [signals, lf_brain, lf_external, freq_dict] = generate_signals(n)
    addpath('ECGSimulation')
    
    signals = cell(n,1);
    
    % Parameters
    nr_sensors  = 123;
    Fs          = 1000; % Hz
    L           = 10000; % Length of signal
    t_trial     = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
    n_trials    = 1;
    t           = (0:(n_trials*length(t_trial)-1)) * (1/Fs);
    n_samples   = length(t);
    i_signal    = 9304; % source in 
    noise_radius  = ["1m env", "2m env", "3m env", "10m env"]; % keys
    freqs       = [12, 16, 28, 40]; % values
    freq_dict = dictionary(noise_radius, freqs);
    env_1m_amp  = 1e-3;
    env_2m_amp  = 1e-2;
    env_3m_amp  = 1e-1;
    env_10m_amp = 1;
    ecg_amp     = 1e-3;
    components  = 3;
    B           = [];
    
    % Generate leadfields
    [lf_brain, n_brain_sources] = generate_leadfield_brain();
    [lf_external, n_external_sources] = generate_leadfield_external();
    
    % For each trail
    for i = 1:n
        %brain signal
        q_signal = brain_signal(t);
        q = zeros(n_brain_sources, n_samples);
        q(i_signal, :) = q_signal;
        B.brain_signal = lf_brain*q;
    
        %sensor noise
        B.sensor_noise = zeros(nr_sensors, n_samples);
        for j = 1:nr_sensors
            [~, ~, sensor_noise_result] = sensor_noise(n_samples);
            B.sensor_noise(j, :) = sensor_noise_result;
        end
    
        %environmental noise
        environmental_noise_1m  = environmental_noise(t, freq_dict("1m env"), env_1m_amp);
        environmental_noise_2m  = environmental_noise(t, freq_dict("2m env"), env_2m_amp);
        environmental_noise_3m  = environmental_noise(t, freq_dict("3m env"), env_3m_amp);
        environmental_noise_10m = environmental_noise(t, freq_dict("10m env"), env_10m_amp);
    
        q = zeros(n_external_sources, n_samples);
        q(2, :) = environmental_noise_1m;
        B.env_1m = lf_external * q;
    
        q = zeros(n_external_sources, n_samples);
        q(3, :) = environmental_noise_2m;
        B.env_2m = lf_external * q;
    
        q = zeros(n_external_sources, n_samples);
        q(4, :) = environmental_noise_3m;
        B.env_3m = lf_external * q;
    
        q = zeros(n_external_sources, n_samples);
        q(5, :) = environmental_noise_10m;
        B.env_10m = lf_external * q;
    
        %biological noise
        [ecg_ft_result, locs] = ecg_ft(components, Fs, L);
        biological_noise = ecg_ft_result * ecg_amp;

        q = zeros(n_external_sources, n_samples);
        q(1, :) = biological_noise;
        B.biological_noise = lf_external * q;
    
        % combined signal
        signal = B.brain_signal + B.sensor_noise + B.biological_noise + B.env_1m + B.env_2m + B.env_3m + B.env_10m;
        signals{i,1} = signal;
    end
    freq_dict("ecg 1") = locs(1,1);
    freq_dict("ecg 2") = locs(1,2);
    freq_dict("ecg 3") = locs(1,3);
    freq_dict("brain signal") = 20;
end