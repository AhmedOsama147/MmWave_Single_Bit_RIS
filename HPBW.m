function theta_HPBW = HPBW(angles,AF)
    [~,peak_idx] = max(AF); theta_peak = angles(peak_idx);
    half_power_dB = -3;
    left_idx = find(AF(1:peak_idx) <= half_power_dB, 1, 'last');
    right_idx = find(AF(peak_idx:end) <= half_power_dB, 1, 'first') + peak_idx - 1;
    
    % Return angles
    left_angle = angles(left_idx);
    right_angle = angles(right_idx);
    if isempty(left_angle) 
        theta_HPBW = 2*abs(right_angle - theta_peak);
    elseif isempty(right_angle)
        theta_HPBW = 2*abs(left_angle - theta_peak);
    else 
        theta_HPBW = abs(left_angle-right_angle);
    end
end