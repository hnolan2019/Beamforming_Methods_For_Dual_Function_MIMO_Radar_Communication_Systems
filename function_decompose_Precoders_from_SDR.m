% ------- Function to decompose the obtained covariance matrices------%

function [Beampattern, success] = function_decompose_Precoders_from_SDR(a, R, h, M, K, Rk)

    % function aims to decompose the obtained covariance matrices from the
    % optimisation problems and determine the precoders
    %
    % Using the precoders the transmit beampattern is calculated
    %
    %
    % Input Arguments 
    % a:              (M x L) - steering vector, 
    % h:              (M x K) - channel vector, 
    % K                       - number of users,
    % M:                      - Number of transmit antennas
    % R              (M x M)  - Overall Covariance Matrix
    % R_k         (M x M x K) - individual covariance matrices for commuication users
    success = false;
    Beampattern = [];

    Wc = function_communication_beamformers(M, K, h, Rk);
    Rc = Wc*Wc';                                % Communication Covariance
    Rr = R-Rc;                                  % Radar Covariance

    try 
       Wr = chol(Rr + eye(M)*1e-9, 'lower');
       W_final = [Wc, Wr];                         % Final Beamforming matrix
       R_final = W_final*W_final';                 % Final Covariance matrix
        
       % ------------- Calculate the SDR optimised beampattern
       Beampattern = real(diag(a'*R_final*a));
       success = true;
    catch
        warning('problem with cholesky');
    end
end