function [ecg_ft_result, locs] = ecg_ft(components, Fs, L)
    T = 1/Fs;           % Sampling period (s)
    t = (0:L-1)*T;      % Time vector
    ecg = ecgSimulation(t);

    ecg_ft = fft(ecg);
    P2 = abs(ecg_ft/L);
    P1 = P2(1:L/2+1);   % Single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1); % Correct amplitude
    
    % Frequency axis
    f_axis = Fs * (0:(L/2)) / L;
    
    % Find Peaks
    [pks, locs] = findpeaks(P1, f_axis, 'MinPeakHeight', 0.1);

    % Plot frequency spectrum with peaks
    % figure;
    % plot(f_axis, P1);
    % xlim([-1 10])
    % hold on;
    % plot(locs, pks, 'ro', 'MarkerFaceColor', 'r'); % Mark peaks
    % title('Magnitude Spectrum with Peaks');
    % xlabel('Frequency (Hz)');
    % ylabel('|X(f)|');
    % grid on;
    % legend('Magnitude Spectrum', 'Peaks');
    
    % Display peak values
    % disp('Peak Frequencies and Magnitudes:');
    % disp(table(locs', pks', 'VariableNames', {'Frequency_Hz', 'Magnitude'}));

    
    phase = angle(ecg_ft);
    phase_deg = rad2deg(phase); % Convert to degrees
    phase_single = phase_deg(1:L/2+1);
    
    ecg_components = zeros(components,3);
    for i = 1:components
        [~, idx] = min(abs(f_axis - locs(1,i)));
        phase_at_target = phase_single(idx);
        ecg_components(i,:) = [locs(1,i), pks(1,i), phase_at_target];
        % fprintf('At %.1f Hz: Magnitude = %.4f, Phase = %.2f degrees\n', locs(1,i), pks(1,i), phase_at_target);
    end
    
    x1 = pks(1,1)*sin(2*pi*locs(1,1)*t+ecg_components(1,3));
    x2 = pks(1,2)*sin(2*pi*locs(1,2)*t+ecg_components(2,3));
    x3 = pks(1,3)*sin(2*pi*locs(1,3)*t+ecg_components(3,3));
    ecg_ft_result = x1+x2+x3;

    % figure;
    % plot(t, ecg_ft_result);
    % title('ECG components');
    % xlabel('Time (s)');
    % ylabel('Amplitude');
    % grid on;
    
end

