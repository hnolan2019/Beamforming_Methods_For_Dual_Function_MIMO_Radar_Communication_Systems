% Simulate the feasbility of the proposed methods at differenct operating points (gamma and K) and export the data for plotting
%
% Feasibility Plot - separated from transmit beampattern as other method
% takes too long
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
ang = -90:0.1:90;                                 % Grid points
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

% ---------------------- Parallel Pool ----------------------%
if isempty(gcp('nocreate'))
    parpool('local', 2);
end

% --------------------- Radar Only --------------------------------------%
des= function_ideal_radar(ang,tgtAz,beamwidth); % Ideal Radar Beampattern 
a= function_steering_vectors(ang, M, d,lambda); % Radar Steering Vector 
A_tgt = a(:, target_cols);                      % Steering vectors to the targets

% Radar Only Optimisation Using MSE beampattern Matching (Used as Benchmark)
radar_only_beampattern_mse = radar_only_optimisation(A_tgt, a, des, L, P, Pt, M, num_pairs, wc);

% --------------------- Communication System ----------------------------%
K_values = [2,4,6];                              % Array of users for communication system
gamma_db_range = 2:2:10;                         % SNIR threshold sweep

%---------------------- Feasibility Plot Memory Allocation --------------%
feasibility_mse = zeros(length(K_values), length(gamma_db_range));
feasibility_final = zeros(length(K_values), length(gamma_db_range));
feasibility_eqr = zeros(length(K_values), length(gamma_db_range));

% ----------------------Equiripple Parameteres -----------------------%
delta_trans = 5;                                % The transition band 
[passband_idx, transition_idx, stopband_idx] = function_define_the_equiripple_regions(ang, beamwidth, delta_trans, tgtAz);

% ------- Determine Feasibility ---------------------------------------%
mc_iterations = 1000;

% Plot the Beampattern for a value of Gamma and different users
for g_idx = 1:length(gamma_db_range)
    gamma_db = gamma_db_range(g_idx);
    gamma = 10^(gamma_db/10);

    for k_idx = 1:length(K_values)
        K = K_values(k_idx);
        theta = linspace(-60, 60, K);


        % Preallocate result arrays before parfor
            success_mse = zeros(1, mc_iterations);
            success_final = zeros(1, mc_iterations);
            success_mxm = zeros(1, mc_iterations);
            
            parfor idx = 1:mc_iterations
                h = function_generate_channel(M, K);
            
                % MSE Optimisation 
                [~, ~, ~, status] = function_SDR_MSE_Joint_Radar_Optimisation(a, A_tgt, des, gamma, h, K, L, P, Pt, M, num_pairs, wc);
                if ~contains(status, 'Solved')
                    continue;
                end
                success_mse(idx) = 1;
            
                % Max-Min with SCA Optimisation 
                [~, ~, ~, ~, status] = function_Max_Min_Successive_Convex_Solution(a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M);
                if ~contains(status, 'Solved')
                    continue;
                end
                success_final(idx) = 1;
            
                % MaxMin Optimisation 
                [~, ~, ~, status] = function_SDR_MaxMin_Optimisation(a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M);
                if ~contains(status, 'Solved')
                    continue;
                end
                success_mxm(idx) = 1;
            end
            
            % Sum after parfor
            successful_runs_mse = sum(success_mse);
            successful_runs_final = sum(success_final);
            successful_runs_mxm = sum(success_mxm);

       
            if (check_runs_were_successful(successful_runs_mse, successful_runs_mxm))
    
                % -------------MSE Average ------------------------------------%
                feasibility_rate_mse = (successful_runs_mse/mc_iterations)*100;
                feasibility_mse(k_idx, g_idx) =(successful_runs_mse /mc_iterations)*100;
    
                % ------------- Final Max-Min Average ------------------------------------%
                feasibility_rate_final_max_min = (successful_runs_final/mc_iterations)*100;
                feasibility_final(k_idx, g_idx)=(successful_runs_final/mc_iterations)*100;
    
                % ---------------Original Max Min Average ---------------------------%
                feasibility_rate_eqr = (successful_runs_mxm/mc_iterations)*100;
                feasibility_eqr(k_idx, g_idx)=(successful_runs_mxm/mc_iterations)*100;
            
            else
                error('All Monte Carlo Iterations failed');
            end
    
        end
    end

export_sidelobe_power(gamma_db_range, K_values, feasibility_mse, feasibility_eqr, feasibility_final);

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

% (2) Export the data to be plotted
function export_sidelobe_power(gamma_db_range, K_values, feasibility_mse, feasibility_eqr, feasibility_final)

    % --- Save data in Data/feasibility_data folder ---
    this_file = mfilename('fullpath');
    this_dir = fileparts(this_file);
    data_dir = fullfile(this_dir, 'Data', 'feasibility_data');

    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    % --- Common first column ---
    Gamma_dB = gamma_db_range(:);

    % --- Build one table per method ---
    T_mse = table(Gamma_dB, 'VariableNames', {'Gamma_dB'});
    T_eqr = table(Gamma_dB, 'VariableNames', {'Gamma_dB'});
    T_final = table(Gamma_dB, 'VariableNames', {'Gamma_dB'});

    for k_idx = 1:length(K_values)
        K = K_values(k_idx);
        col_name = sprintf('K_%d', K);

        T_mse.(col_name) = feasibility_mse(k_idx, :).';
        T_eqr.(col_name) = feasibility_eqr(k_idx, :).';
        T_final.(col_name) = feasibility_final(k_idx, :).';
    end

    % --- Write CSV files ---
    writetable(T_mse, fullfile(data_dir,'feasibility_MSE.csv'));
    writetable(T_eqr, fullfile(data_dir,'feasibility_MaxMin.csv'));
    writetable(T_final, fullfile(data_dir,'feasibility_Proposed_SCA.csv'));
end


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


function successful_runs = check_runs_were_successful(successful_runs_mse, successful_runs_eqr)

    if successful_runs_mse > 0
        successful_runs = true;
    elseif successful_runs_eqr > 0
        successful_runs = true;
    else
        successful_runs = false;
    end
end





