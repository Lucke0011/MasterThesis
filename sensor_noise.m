function [white_noise,pink_noise,sensor_noise_result] = sensor_noise(L)
    % Parameters
    opm_noise_amp = 0.3e-12; % T/Hz^1/2 8e-13 to match real data
    intersect = 0.02; % Pink should intersect white at 10 Hz to match real data
    
    % white_noise
    white_noise_gen = dsp.ColoredNoise(0,'SamplesPerFrame', L);
    white_noise = white_noise_gen()*0.6;
    
    % Pink Noise (1/f noise)
    pink_noise_gen = dsp.ColoredNoise(1.75,'SamplesPerFrame', L);
    pink_noise = pink_noise_gen() * intersect;
    
    sensor_noise_result = (white_noise + pink_noise)*opm_noise_amp;
end

