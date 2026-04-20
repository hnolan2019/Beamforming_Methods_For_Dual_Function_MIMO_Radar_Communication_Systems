function [R, Rk, mydelta, sidelobes_factor_dB_opt, status] =  function_Max_Min_Successive_Convex_Solution(a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M)
    %function to optimise the precoders to match the desired radar beampattern
    % The objective is to design a radar beampattern that maximises the minimum trannsmit power in the passband
    % using sucessive convex approximation

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
    % R        (M x M)         Overall Covariance Matrix
    % R_k      (M x M x K)     individual covariance matrices for commuication users
    % mydelta                  Minimum Transmit Beampattern Power
    % sidelobes_factor_dB_opt  optimal adative sidelobe suppression factor
    % status                   status of the CVX optimisation    
    
    % Fixed iteration parameters
    max_iter = 50;
    tol = 1e-6;
    
    % initial guess for R*
    [R, Rk, mydelta,sidelobes_factor_dB_opt, status] = function_SDR_MaxMin_Optimisation_kappa(a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M);
    R0 = R;
    if ~strcmp(status, 'Solved')
        warning('No initial Guess'); 
        return;
    end
    
    % ------- Calculating delta_0 and epsilon_stop -----%
    A_pass = a(:,passband_idx);
    A_stop = a(:,stopband_idx);
    
    % ------- Calculating the stopband and passband ------- %
    P_pass0 = real(diag(A_pass'*R0*A_pass));        
    P_stop0 = real(diag(A_stop'*R0*A_stop));        
    

    % -------- Start Iteration 1 -----------------------------%
    delta_n = min(P_pass0);                         % delta^(r)
    eps_n = max(P_stop0);                           % epsilon_stop^(r)
    eps_n = max(eps_n, 1e-6);                       % avoid divide by zero
    obj_prev = delta_n/eps_n;
    
    %-----------------SCA iterations ------------------------$
    for iter = 1:max_iter
        
        % partial derivatives each iteration
        grad_h_delta = (delta_n + 1)/(2*eps_n);
        grad_h_eps = -(delta_n + 1)^2/(4*eps_n^2);
        h_const = (delta_n + 1)^2/(4 *eps_n);
        
        cvx_clear;
        cvx_begin sdp quiet
            variable R(M,M) hermitian semidefinite
            variable Rk(M,M,K) hermitian semidefinite
            variable mydelta(1) nonnegative
            variable epsilon_stop(1) nonnegative
            
            P_pass = real(diag(A_pass'*R*A_pass));
            P_stop = real(diag(A_stop'*R*A_stop));
            
            g = 0.25*quad_over_lin(mydelta - 1, epsilon_stop);
            
            % first order taylor expansion 
            h_tilda = h_const + grad_h_delta*(mydelta-delta_n) + grad_h_eps*(epsilon_stop-eps_n);
            
            % optimisation problem 
            minimize(g-h_tilda);
              
            % ------ Constraints --------%
            subject to
                epsilon_stop >= 1e-6;
                P_pass >= mydelta;
                P_stop <= epsilon_stop;
                
                % ----Cross Correlation ------%
                for p = 1:P-1
                    for q = p+1:P
                        abs(A_tgt(:,q)'*R*A_tgt(:,p)) <= 0.05;
                    end
                end
                
                diag(R) == Pt / M;
                
                % --- SINR constraints ----%
                for k = 1:K
                    hk = h(:,k);
                    signal_power = real(hk'*Rk(:,:,k)*hk);
                    interference_power = real(hk'*R*hk) - signal_power;
                    signal_power >= gamma*(interference_power + 1);
                end
                R - sum(Rk, 3) == semidefinite(M);
        cvx_end
        
        % Check if CVX failed
        if ~strcmp(cvx_status, 'Solved')
            warning('SCA iteration %d failed: %s', iter, cvx_status);
            break;
        end
        
        % --- Update operating point ---
        delta_n = mydelta;
        eps_n = max(epsilon_stop, 1e-6);
        obj_curr = delta_n / eps_n;
        
        
        % --- Check convergence ---
        if abs(obj_curr-obj_prev) / (1 + abs(obj_prev)) < tol
            fprintf('Converged at iteration %d\n', iter);
            break;
        end
        obj_prev = obj_curr;
    end
    
    % compute the optimal suppression value
    kappa_opt = mydelta/epsilon_stop;
    sidelobes_factor_dB_opt = 10*log10(kappa_opt);
    status = cvx_status;
    fprintf('ratio delta/eps:  %.6e\n', mydelta/epsilon_stop);
end