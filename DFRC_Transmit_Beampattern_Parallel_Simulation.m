% --- numerically evaluate the MIMO radar transmit beam patterns -------%
% Simulate 
% Paper MSE
% Conference Paper Solution
% Final Solution
%
%
% -------------------Fixed Constsants ----------------------------------%
% (1) Tx antennas are equally spaced: d = half wavelength
% (2) The Number of Transmit Antennas (M): 16
% (3) Total Transmit Power (Pt): 16dB
% (4) 1 Target directions from basestation: 0 degrees
% (5) Ideal beamwidth: 10 degrees
% (6) Grid points obtained by unform sampling the range [-90 to 90] with
%     resolution 0.1

clear all; close all; clc

% ---------------------- Fixed Constants--------------------%
fc = 77e9; c = physconst('LightSpeed'); lambda = c/fc;
d = lambda/2;                                   % Antenna Spacing
ang = -90:0/1:90;                               % Grid points
M=16;                                           % Number of transmit Antennas
Pt = M;                                         % total tarnsmit power ( per antenna power = 1)
tgtAz = [-40, 0, 40];                           % Location of radar targets
P = length(tgtAz);                              % The number of targets
L = length(ang);                                % The number of grid points
tgtRng = [5.31e3, 6.23e3, 5.7e3];               % Ranges of the targets
beamwidth = 10;                                 % fixed beam width
num_pairs = (P^2 - P) / 2;                      % Number of unique radar target pairs
target_cols = targetsCol(P, ang, tgtAz);        % index of the targets in ang
wc = 1;                                         % Weighting factor

% Set the parallel code
if isempty(gcp('nocreate'))
    parpool('local');
end


% --------------------- Radar Only --------------------------------------%
des= function_ideal_radar(ang,tgtAz,beamwidth); % Ideal Radar Beampattern 
a= function_steering_vectors(ang, M, d,lambda); % Radar Steering Vector 
A_tgt = a(:, target_cols);                      % Steering vectors to the targets

% Radar Only Optimisation Using MSE beampattern Matching (Used as Benchmark)
radar_only_beampattern_mse = radar_only_optimisation(A_tgt, a, des, L, P, Pt, M, num_pairs, wc);

% --------------------- Communication System ----------------------------%
K_values = [2,4,6];                              % Array of users for communication system
gamma_db_range = 10;                         % SNIR threshold sweep



% ----------------------Equiripple Parameteres -----------------------%
delta_trans = 5;                                % The transition band 
[passband_idx, transition_idx, stopband_idx] = function_define_the_equiripple_regions(ang, beamwidth, delta_trans, tgtAz);

% --------------------- Joint Radar-Communication Optimisation
mc_iterations = 1000;

% Plot the Beampattern for a value of Gamma and different users
for g_idx = 1:length(gamma_db_range)
    gamma_db = gamma_db_range(g_idx);
    gamma = 10^(gamma_db/10);

    fprintf("The Gamma Value is: %f", gamma);


    for k_idx = 1:length(K_values)
        K = K_values(k_idx);
        theta = linspace(-60, 60, K);

        fprintf("The Number of users is: %d", K);

        % ----- Allocate Memory for each DFRC Optimisation Approach ---------%
        MSE_beampattern_sum = zeros(L, 1);
        finalMax_min_beampattern_sum = zeros(L, 1);
        max_min_beampattern_sum = zeros(L,1);

        alpha_opt_sum = 0;

        % Preallocate per-iteration outputs
        mse_bp_all = zeros(L, mc_iterations);
        mxm_sca_bp_all = zeros(L, mc_iterations);
        mxm_bp_all = zeros(L, mc_iterations);
        
        alpha_all = zeros(1, mc_iterations);
        
        success_mse = false(1, mc_iterations);
        success_final = false(1, mc_iterations);
        success_mxm = false(1, mc_iterations);

            parfor idx = 1:mc_iterations
                local_mse_bp = zeros(L,1);
                local_final_bp = zeros(L,1);
                local_mxm_bp = zeros(L,1);
            
                local_alpha = 0;
            
                local_success_mse = false;
                local_success_final = false;
                local_success_mxm = false;
        
            
                % ----------generate the channel ----------%
                h = function_generate_channel(M,K);
            
                % ---- Simulation MSE Beampattern Matching ----%
                [R, Rk, alpha_opt, status] = function_SDR_MSE_Joint_Radar_Optimisation(a, A_tgt, des, gamma, h, K, L, P, Pt, M, num_pairs, wc);
            
                % -------- If MSE solved correctly proceed with decomposition ----%
                if contains(status, 'Solved')
                    [MSE_beampattern, success] = function_decompose_Precoders_from_SDR(a, R, h, M, K, Rk);
                    if success
                        local_success_mse = true;
                        local_alpha = alpha_opt;
                        local_mse_bp = MSE_beampattern(:);
                    end
                end
            
                % ---- Simulation Max-Min Successive Convex Solution ----
                [R, Rk, mydelta, sidelobes_factor_dB_opt, status] = function_Max_Min_Successive_Convex_Solution(a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M);
            
                 % -------- If Max-Min solved correctly proceed with decomposition ----%
                if contains(status, 'Solved')
                    [Final_Max_Min_beampattern, success] = function_decompose_Precoders_from_SDR(a, R, h, M, K, Rk);
                    if success
                        local_success_final = true;
                        local_final_bp = Final_Max_Min_beampattern(:);
                    end
                end
            
                % ---- Simulation original Max Min Beampattern Matching ----
                [R, Rk, delta, status] = function_SDR_MaxMin_Optimisation( a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M);
            
                 % -------- If Max-Min solved correctly proceed with decomposition ----%
                if contains(status, 'Solved')
                    [mxm_beampattern, success] = function_decompose_Precoders_from_SDR(a, R, h, M, K, Rk);
                    if success
                        local_success_mxm = true;
                        local_mxm_bp = mxm_beampattern(:);
                    end
                end
            
                mse_bp_all(:,idx) = local_mse_bp;
                mxm_sca_bp_all(:,idx) = local_final_bp;
                mxm_bp_all(:,idx) = local_mxm_bp;
            
                alpha_all(idx) = local_alpha;
            
                success_mse(idx) = local_success_mse;
                success_final(idx) = local_success_final;
                success_mxm(idx) = local_success_mxm;
            end

        successful_runs_mse = sum(success_mse);
        successful_runs_final = sum(success_final);
        successful_runs_mxm = sum(success_mxm);
        
        if check_runs_were_successful(successful_runs_mse, successful_runs_mxm)
        
            feasibility_rate_mse = 100 * successful_runs_mse/mc_iterations;
            fprintf('\n Simulation feasibility rate: %f', feasibility_rate_mse);
        
            MSE_Beampattern_Avg = sum(mse_bp_all(:, success_mse), 2)/successful_runs_mse;
            alpha_opt = sum(alpha_all(success_mse))/successful_runs_mse;
        
            feasibility_rate_sca_max_min = 100*successful_runs_final/mc_iterations;
            fprintf('\n Simulation feasibility rate: %f', feasibility_rate_sca_max_min);
            SCA_MXM_Beampattern_Avg = sum(mxm_sca_bp_all(:, success_final), 2)/successful_runs_final;
        
            feasibility_rate_mxm = 100*successful_runs_mxm/mc_iterations;
            fprintf('\n Simulation feasibility rate: %f', feasibility_rate_mxm);
            MAX_Beampattern_Avg = sum(mxm_bp_all(:, success_mxm), 2)/successful_runs_mxm;
        
            export_radar_beampattern(MSE_Beampattern_Avg, SCA_MXM_Beampattern_Avg,MAX_Beampattern_Avg, ang, radar_only_beampattern_mse, K, gamma, des, alpha_opt);
        else
            error('All Monte Carlo Iterations failed');
        end
    end
end



% ------- FUNCTIONS -----------------------------------%

% (1) Radar Only Optimisation using MSE Beampattern matching
function radar_only_beampattern = radar_only_optimisation(A_tgt, a, des, L, P, Pt, M, num_pairs, wc)
    cvx_clear;
    cvx_begin sdp quiet
        variable R_r(M, M) hermitian semidefinite 
        variable alpha_opt
    
        P_pattern_opt = real(diag(a'*R_r*a));    
        loss_beampattern_error = (1/L)*sum_square(alpha_opt* des' - P_pattern_opt);
    
        loss_crosscorr = 0;
        for p = 1:P-1
            for q = p+1:P
                loss_crosscorr = loss_crosscorr + sum_square_abs(A_tgt(:,q)'*R_r*A_tgt(:,p));
            end
        end
        loss_crosscorr = (1/num_pairs)*loss_crosscorr;
        minimize(loss_beampattern_error + wc*loss_crosscorr)
    
        % constraints
        subject to
            diag(R_r) == Pt/M;
    
    cvx_end

    radar_only_beampattern = real(diag(a'*R_r*a));
end

% --------------- Function determines the index of the target array ------%
function target_cols = targetsCol(P, ang, tgtAz)
 target_cols = zeros(size(P));            % Array for index of the targets in (ang)
     for i = 1:P
        idx = find(ang == tgtAz(i), 1);
        if isempty(idx)
            error('TARGETS MUST BE KNOWN');
        end
        target_cols(i) = idx;       
     end
end


function export_radar_beampattern(MSE_Beampattern_Avg, SCA_MXM_Beampattern_Avg, MAX_Beampattern_Avg, ang, radar_only_beampattern_mse, K, gamma, des, alpha_opt)
    gamma_db = 10*log10(gamma);
    
    % ---Sava data in the Data/Transmit_Beampattern folder for plotting
    this_file = mfilename('fullpath');
    this_dir = fileparts(this_file);
    data_dir = fullfile(this_dir, 'Data/Transmit_Beampattern');

    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    % --- Force column vectors so table
    ang  = ang(:);
    MSE_col  = MSE_Beampattern_Avg(:);
    FinalMXM_col  = SCA_MXM_Beampattern_Avg(:);
    MaxMin_col = MAX_Beampattern_Avg(:);
    Ideal_col = alpha_opt*des(:);
    RadarOnly_col = radar_only_beampattern_mse(:);

    % --- Build the table 
    T = table(ang, MSE_col, FinalMXM_col, MaxMin_col, Ideal_col, RadarOnly_col, ...
        'VariableNames', { ...
            'Azimuth_deg', ...
            'MSE', ...
            'FinalMaxMin', ...
            'MaxMin', ...
            'Ideal', ...
            'RadarOnly'});

    % --- Filename encodes the operating point ---
    fname = sprintf('beampattern_results_%d_users_%ddB.csv', K, round(gamma_db));
    fpath = fullfile(data_dir, fname);

    writetable(T, fpath);
end



function successful_runs = check_runs_were_successful(successful_runs_mse, successful_runs_eqr)

    if successful_runs_mse > 0
        successful_runs = true;
    elseif successful_runs_eqr > 0
        successful_runs = true;
    else
        successful_runs = false;
    end
end





