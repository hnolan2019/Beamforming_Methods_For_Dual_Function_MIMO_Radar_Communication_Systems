% ----------------- Function to perform the CVX Optimisation --------------%
function [R, Rk, alpha_opt, status] = function_SDR_MSE_Joint_Radar_Optimisation(a, A_tgt, des, gamma, h, K, L, P, Pt, M,num_pairs, wc)
    %function to optimise the precoders to match the desired radar beampattern
    % The objective is to design a radar beampattern that minimises the
    % mean squared error between the obtained radar beampattern and the
    % ideal radar beampattern

    % Input Arguments 
    % a:             (M x L) - steering vector, 
    % A_tgt:         (M x P) - index of the steering vectors in ang,
    % des            (1 x L) - The ideal radar beampattern
    % gamma:                 - Minimum required SINR at each user
    % h:             (M x K) - channel vector, K = number of users,
    % K                        The number of users in the system
    % L                        The number of grid points in ang
    % P:                       Number of radar targets
    % Pt:                      Total power transmitted from the basestation
    % M:                       Number of transmit antennas
    % num_oairs:               Number of unique pairs of radar targets
    % wc:                      Weighting factor to decide on importance 
    %                          between cross-correlation and beampattern error minimisation
    %
    % Outputs
    % R        (M x M)      Overall Covariance Matrix
    % R_k      (M x M x K)  individual covariance matrices for commuication users
    % alpha_opt             optimial value for desired beampattern
    % status                status of the CVX optimisation
  
  
  
% --------------Joint Radar - Communication Optimisation Problem -------%
    cvx_clear;
    cvx_begin sdp quiet
        variable R(M, M) hermitian semidefinite         % total covariance
        variable Rk(M,M,K) hermitian semidefinite       % per-user covariance matrix
        variable alpha_opt(1) nonnegative
    
        % --------------------- Beampattern with MSE ----------------------%
        P_pattern_opt = real(diag(a'*R*a));             % Radar Pattern
        loss_beampattern_error = (1/L)*sum_square(alpha_opt* des' - P_pattern_opt);
    
        % --------------------- Cross Correlation Function ----------------%
        loss_crosscorr = 0;
        for p = 1:P-1
            for q = p+1:P
                loss_crosscorr = loss_crosscorr + sum_square_abs(A_tgt(:,q)'*R*A_tgt(:,p));
            end
        end
        loss_crosscorr = (1/num_pairs)*loss_crosscorr;
    
        % --------------------- Optimisation Function ---------------------%
        minimize(loss_beampattern_error + wc * loss_crosscorr);
    
        % ---------------------  Optimisation Constraints -----------------%
        subject to
            diag(R) == Pt/M;
    
            for k = 1:K
                hk = h(:,k);
                signal_power = real(hk'*Rk(:,:,k)*hk);
                interference_power = real(hk'*R*hk) - signal_power;
                signal_power >= gamma * (interference_power +1);
            end 
            sum(Rk, 3) <= R;
    cvx_end
    status = cvx_status;
end
