function dental_noise_result = dental_noise(t, f, amp)
    dental_noise_result = amp .* sin(2*pi*f*t);
end

