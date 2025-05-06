function signals = empty_room_signals(n)

    signals = cell(n,1);
    
    % Parameters
    nr_sensors  = 246;
    Fs          = 1000; % Hz
    t_trial     = -0.5:(1/Fs):(9.5-1/Fs); % -0.5:+9.499 s
    n_trials    = 1;
    t           = (0:(n_trials*length(t_trial)-1)) * (1/Fs);
    n_samples   = length(t);
    noise_radius  = ["1m", "2m", "3m", "10m"]; % keys
    freqs       = [12, 16, 28, 40]; % values
    noise_freq_dict = dictionary(noise_radius, freqs);
    env_1m_amp  = 1e-2;
    env_2m_amp  = 1e-1;
    env_3m_amp  = 1;
    env_10m_amp = 1e2;
    B           = [];
    
    % Generate leadfields
    [lf_external, n_external_sources] = generate_leadfield_external();
    
    % For each trail
    for i = 1:n
    
        %sensor noise
        B.sensor_noise = zeros(nr_sensors, n_samples);
        for j = 1:nr_sensors
            [~, ~, sensor_noise_result] = sensor_noise(n_samples);
            B.sensor_noise(j, :) = sensor_noise_result;
        end
    
        %environmental noise
        environmental_noise_1m  = environmental_noise(t, noise_freq_dict("1m"), env_1m_amp);
        environmental_noise_2m  = environmental_noise(t, noise_freq_dict("2m"), env_2m_amp);
        environmental_noise_3m  = environmental_noise(t, noise_freq_dict("3m"), env_3m_amp);
        environmental_noise_10m = environmental_noise(t, noise_freq_dict("10m"), env_10m_amp);
    
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
        signal = B.sensor_noise + B.env_1m + B.env_2m + B.env_3m + B.env_10m;
        signals{i,1} = signal;
    end
end

