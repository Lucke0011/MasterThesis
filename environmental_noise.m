function environmental_noise_result = environmental_noise(t, f, amp)
    phase_offset = 2 * pi * rand(); % simulate variationsin the inter trial interval 
    environmental_noise_result = amp .* sin(2*pi*f*t + phase_offset); %A*m
end

