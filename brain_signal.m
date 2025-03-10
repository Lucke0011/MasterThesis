function q_signal = brain_signal(t)
    Q_signal = 50e-9; %A*m
    f_signal = 20; % Hz
    q_signal = Q_signal*cos(2*pi*f_signal*t);
end

