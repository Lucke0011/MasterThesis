function [signals, lf_brain, lf_external, freq_dict] = signals_no_ecg()
    
    % Parameters
    nr_sensors  = 123;
    Fs          = 1000; % Hz
    L           = 10000; % Length of signal
    t_trial     = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
    n_trials    = 1;
    t           = (0:(n_trials*length(t_trial)-1)) * (1/Fs);
    n_samples   = length(t);
    sources = [15240, 13921, 11965, 9101, 7412, 6109, 4782, 3899, 2294, 1001];
    signals = cell(length(sources),1);
    noise_radius  = ["1m env", "2m env", "3m env", "10m env"]; % keys
    freqs       = [12, 16, 28, 40]; % values
    freq_dict = dictionary(noise_radius, freqs);
    env_1m_amp  = 1e-3;
    env_2m_amp  = 1e-2;
    env_3m_amp  = 1e-1;
    env_10m_amp = 1;
    B           = [];
    
    % Generate leadfields
    [lf_brain, n_brain_sources] = generate_leadfield_brain();
    [lf_external, n_external_sources] = generate_leadfield_external();
    
    % For each source
    for i = 1:length(sources)
        %brain signal
        q_signal = brain_signal(t);
        q = zeros(n_brain_sources, n_samples);
        q(sources(i), :) = q_signal;
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
    
        % combined signal
        signal = B.brain_signal + B.sensor_noise + B.env_1m + B.env_2m + B.env_3m + B.env_10m;
        signals{i,1} = signal;
    end
    freq_dict("brain_signal") = 20;
end
