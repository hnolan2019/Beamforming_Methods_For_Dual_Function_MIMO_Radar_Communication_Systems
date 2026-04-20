% ------- Function to generate the radar steering vectors ------%
function a = function_steering_vectors(ang, M, d, lambda)

% Inputs - grid points, number of transmit antennas, antenna spacing (d), 
%          radar waveform
% Outputs - radar steering vector

    m  = (0:M-1)';
    theta_rad = deg2rad(ang);
    a  = exp(-1j*(2*pi)*d/lambda*m*sin(theta_rad));
end