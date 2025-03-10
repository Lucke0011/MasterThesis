function psd_diff_result = psd_diff(signal_before, signal_after, noise_freq_dict, Fs)
    psd_diff_result = cell(length(keys(noise_freq_dict)), 1);
    
    % Before HFC
    freqs = keys(noise_freq_dict);
    
    for i = 1:length(freqs)
        freq = noise_freq_dict(freqs(i));
    
        [pxx_signal, f] = pwelch(signal_before', [], [], 0:0.2:500, Fs); % 0.2 Hz
        
        % Highest value measured
        psd_before = 0;
        for j = 1:123
            temp = interp1(f, pxx_signal(:,j), freq);
            if temp > psd_before
                psd_before = temp;
            end
        end
        %db_before = 10*log10(psd_before);
        fprintf('The dB before value of %s at f = %.2f Hz is %f dB \n', freqs(i), freq, psd_before);
        
        % After HFC
        [pxx_signal, f] = pwelch(signal_after', [], [], 0:0.2:500, Fs);
        
        % Highest value measured
        psd_after = 0;
        for j = 1:123
            temp = interp1(f, pxx_signal(:,j), freq);
            if temp > psd_after
                psd_after = temp;
            end
        end
        %db_after = 10*log10(psd_after);
        fprintf('The dB after value of %s at f = %.2f Hz is %f dB \n', freqs(i), freq, psd_after);
        
        % Noise difference
        hfc_psd_diff = psd_before / psd_after; % In psd use /
        fprintf('The dB difference of %s at f = %.2f Hz is %f dB \n', freqs(i), freq, hfc_psd_diff);

        psd_diff_result{i} = [freqs(i), freq, psd_before, psd_after, hfc_psd_diff];
    end
end

