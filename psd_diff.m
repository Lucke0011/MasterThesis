function [max_diff, max_one_channel] = psd_diff(signal_before, data_after, freq_dict, Fs, n)
    % max_diff = cell(length(keys(freq_dict)), 1);
    % max_one_channel = cell(length(keys(freq_dict)), 1);

    max_diff = zeros(length(keys(freq_dict)), 1);
    max_one_channel = zeros(length(keys(freq_dict)), 1);


    % Before
    freqs = keys(freq_dict);
    for i = 1:length(freqs)
        freq = freq_dict(freqs(i));
    
        [pxx_signal, f] = pwelch(signal_before', [], [], 0:0.2:500, Fs); % 0.2 Hz
        
        % Highest value measured
        psd_before = 0;
        index = 0;
        for j = 1:123
            temp = interp1(f, pxx_signal(:,j), freq);
            if temp > psd_before
                psd_before = temp;
                index = j;
            end
        end
        db_before = 10*log10(psd_before);
        % fprintf('The dB before value of %s at f = %.2f Hz is %f dB \n', freqs(i), freq, db_before);
            
        for k = 1:n
            % After
            [pxx_signal, f] = pwelch(data_after{k}', [], [], 0:0.2:500, Fs);
            
            % Highest value measured
            psd_after = 0;
            for j = 1:123
                temp = interp1(f, pxx_signal(:,j), freq);
                if temp > psd_after
                    psd_after = temp;
                end
            end
            db_after = 10*log10(psd_after);
            % fprintf('The dB after value of %s at f = %.2f Hz is %f dB \n', freqs(i), freq, db_after);
            
            % PSD/dB difference
            psd_diff = psd_before / psd_after;
            db_diff = db_before - db_after;
            % fprintf('The dB difference of %s at f = %.2f Hz is %f dB \n', freqs(i), freq, db_diff);
    
            % PSD/dB difference one channel
            psd_diff_one_channel = psd_before / interp1(f, pxx_signal(:,index), freq);
            db_diff_one_channel = db_before - (10*log10(interp1(f, pxx_signal(:,index), freq)));
            % fprintf('The dB difference of %s at f = %.2f Hz is %f dB \n', freqs(i), freq, db_diff_one_channel);
    
            % max_diff(i, k)        = psd_diff;
            % max_one_channel(i, k) = psd_diff_one_channel;
            max_diff(i, k)        = db_diff;
            max_one_channel(i, k) = db_diff_one_channel;
        end
    end
end

