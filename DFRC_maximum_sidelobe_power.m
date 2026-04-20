% --- Compare the Maximum Sidelobe Power Vs SINR  -------%
% 
% - Expect the MSE to have a greater sidelobe power relative to the Max-Min
% optimisation problems
%
% -------------------Fixed Constsants ----------------------------------%
% (1) Tx antennas are equally spaced: d = half wavelength
% (2) The Number of Transmit Antennas (M): 16
% (3) Total Transmit Power (Pt): 16dB
% (4) 1 Target directions from basestation: 0 degrees
% (5) Ideal beamwidth: 10 degrees
% (6) Grid points obtained by unform sampling the range [-90 to 90] with
%     resolution 0.5

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

% --------------------- Radar Only --------------------------------------%
des= function_ideal_radar(ang,tgtAz,beamwidth); % Ideal Radar Beampattern 
a= function_steering_vectors(ang, M, d,lambda); % Radar Steering Vector 
A_tgt = a(:, target_cols);                      % Steering vectors to the targets

% Radar Only Optimisation Using MSE beampattern Matching (Used as Benchmark)
radar_only_beampattern_mse = radar_only_optimisation(A_tgt, a, des, L, P, Pt, M, num_pairs, wc);

% --------------------- Communication System ----------------------------%
K_values = [2,4,6];
gamma_db_range = 2:2:10;
mc_iterations = 1000;

% ---------------------- Equiripple Regions -----------------%
delta_trans = 5;
[passband_idx, ~, stopband_idx] = function_define_the_equiripple_regions(...
    ang, beamwidth, delta_trans, tgtAz);

% ---------------------- Parallel Pool ----------------------%
if isempty(gcp('nocreate'))
    parpool('local', 2);
end

% ---------------------- Results Storage --------------------%
MSE_Radar_Performance = NaN(length(K_values), length(gamma_db_range));
MXM_Radar_Performance = NaN(length(K_values), length(gamma_db_range));
MXM_SCA_Radar_Performance = NaN(length(K_values), length(gamma_db_range));

% --- Sidelobe Power Storage ---%
MSE_Sidelobe_Power = NaN(length(K_values), length(gamma_db_range));
MXM_Sidelobe_Power = NaN(length(K_values), length(gamma_db_range));
MXM_SCA_Sidelobe_Power = NaN(length(K_values), length(gamma_db_range));

for k_user = 1:length(K_values)
    K = K_values(k_user);
    theta = linspace(-60, 60, K);

    if ~validate_parameters(M, K, theta, tgtAz)
        error('System parameters are not valid');
    end

    for g_idx = 1:length(gamma_db_range)
        gamma_db = gamma_db_range(g_idx);
        gamma = 10^(gamma_db/10);

        fprintf('\nK=%d | gamma=%.0f dB', K, gamma_db);

        % Preallocate per-iteration accumulators
        success_mse = zeros(1, mc_iterations);
        success_mxm = zeros(1, mc_iterations);
        success_mxm_sca = zeros(1, mc_iterations);

        % --- Sidelobe power accumulators ---%
        mse_sidelobe_acc = zeros(1, mc_iterations);
        mxm_sidelobe_acc = zeros(1, mc_iterations);
        mxm_sca_sidelobe_acc = zeros(1, mc_iterations);

        parfor idx = 1:mc_iterations
            h = function_generate_channel(M, K);

            % --- MSE ---%
            [R, Rk, ~, status] = function_SDR_MSE_Joint_Radar_Optimisation(a, A_tgt, des, gamma, h, K, L, P, Pt, M, num_pairs, wc);
            if ~contains(status, 'Solved'); continue; end

            [bp, ok] = function_decompose_Precoders_from_SDR(a, R, h, M, K, Rk);
            if ~ok; continue; end
            mse_sidelobe_acc(idx) = mean(bp(stopband_idx));
            success_mse(idx) = 1;

            % --- Max-Min  ---%
            [R, Rk, ~, status] = function_SDR_MaxMin_Optimisation(...
                a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M);
            if contains(status, 'Solved')
                [bp, ok] = function_decompose_Precoders_from_SDR(a, R, h, M, K, Rk);
                if ok
                    mxm_sidelobe_acc(idx) = mean(bp(stopband_idx));
                    success_mxm(idx) = 1;
                end
            end

            % --- Max-Min SCA  ---%
            [R, Rk, ~, ~, status] = function_Max_Min_Successive_Convex_Solution(...
                a, A_tgt, h, K, gamma, passband_idx, stopband_idx, P, Pt, M);
            if contains(status, 'Solved')
                [bp, ok] = function_decompose_Precoders_from_SDR(a, R, h, M, K, Rk);
                if ok
                    mxm_sca_sidelobe_acc(idx) = mean(bp(stopband_idx));
                    success_mxm_sca(idx) = 1;
                end
            end
        end

        % --- Average results ---%
        n_mse = sum(success_mse);
        n_mxm = sum(success_mxm);
        n_mxm_sca = sum(success_mxm_sca);

        if n_mse > 0
            MSE_Sidelobe_Power(k_user, g_idx) = 10*log10(sum(mse_sidelobe_acc) / n_mse);
            fprintf('  MSE: %.0f%%', (n_mse/mc_iterations)*100);
        end

        if n_mxm > 0
            MXM_Sidelobe_Power(k_user, g_idx) = 10*log10(sum(mxm_sidelobe_acc) / n_mxm);
            fprintf('  Max-Min: %.0f%%', (n_mxm/mc_iterations)*100);
        end

        if n_mxm_sca > 0
            MXM_SCA_Sidelobe_Power(k_user, g_idx) = 10*log10(sum(mxm_sca_sidelobe_acc) / n_mxm_sca);
            fprintf('  SCA: %.0f%%', (n_mxm_sca/mc_iterations)*100);
        end
    end
end

% ---------------------- Plot --------------------------------%
export_sidelobe_power(gamma_db_range, K_values, MSE_Sidelobe_Power, MXM_Sidelobe_Power, MXM_SCA_Sidelobe_Power)

% ---------- Functions -----------------------%
function isValid = validate_parameters(M, K, theta, tgtAz)
    arguments
        M (1,1) {mustBeInteger, mustBePositive} 
        K (1,1) {mustBeInteger, mustBePositive}
        theta (1,:)  {mustBeReal}
        tgtAz (1,:) {mustBeReal}
    end

    isValid = true;
    if length(theta) ~= K
        isValid = false;
    end
    if any(abs(theta) > 90)
        isValid = false;
    end
    if any(abs(tgtAz) > 90)
        isValid = false;
    end
end


function export_sidelobe_power(gamma_db_range, K_values, MSE_pow, MXM_pow, SCA_pow)
 % ---Sava data in the Data/Transmit_Beampattern folder for plotting
    this_file = mfilename('fullpath');
    this_dir = fileparts(this_file);
    data_dir = fullfile(this_dir, 'Data/Sidelobes');

    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    rows = [];
    for k_idx = 1:length(K_values)
        K = K_values(k_idx);
        for g_idx = 1:length(gamma_db_range)
            gamma_db = gamma_db_range(g_idx);

            rows = [rows; 
                K, gamma_db, MSE_pow(k_idx, g_idx), MXM_pow(k_idx, g_idx), SCA_pow(k_idx, g_idx)];
        end
    end
     T = array2table(rows, ...
        'VariableNames', {'K', 'Gamma_dB', 'MSE_Sidelobe_dB', 'MaxMin_Sidelobe_dB', 'SCA_Sidelobe_dB'});

    % --- Filename encodes the operating point ---
    fname = sprintf('sidelobe_power_data.csv', K, round(gamma_db));
    fpath = fullfile(data_dir, fname);

    writetable(T, fpath);
end


function target_cols = targetsCol(P, ang, tgtAz)
    target_cols = zeros(size(P));
    for i = 1:P
        idx = find(ang == tgtAz(i), 1);
        if isempty(idx)
            error('TARGETS MUST BE KNOWN');
        end
        target_cols(i) = idx;
    end
end

% --------- CVX Optimisation for Radar Only --------------- %
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
        minimize(loss_beampattern_error + wc * loss_crosscorr)
    
        subject to
            diag(R_r) == Pt/M;
    
    cvx_end

    radar_only_beampattern = real(diag(a'*R_r*a));
end