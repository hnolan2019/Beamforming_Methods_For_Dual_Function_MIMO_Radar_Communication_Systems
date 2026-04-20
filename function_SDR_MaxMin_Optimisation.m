% ----- MaxMin  Beampattern Matching Optimisation with fixed sidelobe suppression --------------%
function [R, Rk, mydelta, status] = function_SDR_MaxMin_Optimisation(a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M)

    % function to optimise the precoders to match the desired radar beampattern
    % The objective is to design a radar beampattern that maximises the minimum trannsmit power in the passband

    % Input Arguments 
    % a:             (M x L) - steering vector, 
    % A_tgt:         (M x P) - index of the steering vectors in ang,
    % h:             (M x K) - channel vector, K = number of users,
    % gamma:                 - Minimum required SINR at each user
    % passband_idx:  (1 x L) - index of the passbands in ang, 
    % stopoband idx: (1 x L) - index of the stopbands in ang
    % P:                       Number of radar targets
    % Pt:                      Total power transmitted from the basestation
    % M:                       Number of transmit antennas
    %
    % Outputs
    % R        (M x M)      Overall Covariance Matrix
    % R_k      (M x M x K)  individual covariance matrices for commuication users
    % mydelta                 Minimum Transmit Beampattern Power
    % status                status of the CVX optimisation
    
    A_pass = a(:, passband_idx);
    A_stop = a(:, stopband_idx);
    sidelobes_factor_dB = 12;
    kappa = 10^(-sidelobes_factor_dB/10);
    epsillon = 0.05;
   
    % ------ Clear CVX before optimisation ----------%
    cvx_clear;

    % ----------  Optimisation ------------%
    cvx_begin sdp quiet
        variable R(M,M) hermitian semidefinite
        variable Rk(M,M,K) hermitian semidefinite
        variable mydelta(1) nonnegative  % Variable to maximise

        %---------- Calculate power in passband and stopband----------%
        P_pass = real(diag(A_pass'*R*A_pass)); 
        P_stop = real(diag(A_stop'*R*A_stop)); 

        maximize(mydelta)
    
        subject to
       % ----------------- Radar Constraints ------------%

       P_pass >= mydelta;
       P_stop <= kappa*mydelta;
        % --------------------- Cross Correlation Function ----------------%
        for p = 1:P-1
            for q = p+1:P
                abs(A_tgt(:,q)'*R*A_tgt(:,p)) <= epsillon; 
            end
        end

        diag(R) == Pt / M;              

        % ------------communication constraint ---------------------%
             for k = 1:K                             
                 hk = h(:,k);
                 signal_power = real(hk'*Rk(:,:,k)*hk);
                 interference_power = real(hk'*R*hk)-signal_power;
                 signal_power >= gamma * (interference_power +1);
             end

            R - sum(Rk, 3) == semidefinite(M);
    cvx_end

    status = cvx_status;
end
