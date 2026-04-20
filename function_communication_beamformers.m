% ---------Calculating the Communication Beamformers
function Wc = function_communication_beamformers(M, K, h, Rk)
    % function aims to extract the communication precoders
    %
    % Input Arguments 
    % M:                      - Number of transmit antennas
    % K                       - number of users,
    % M:                      - Number of transmit antennas
    % h:              (M x K) - channel vector,
    % R_k         (M x M x K) - individual covariance matrices for commuication users
    %
    Wc = zeros(M, K);
    
    for k = 1:K
        hk = h(:, k);
        Rk_for_k = Rk(:, :, k);
        wk = ((hk'*Rk_for_k*hk)^(-0.5))*(Rk_for_k*hk);
        Wc(:, k) = wk;
    end
end
