clear; close all; clc;

% ------ Fixed Constants
c = physconst('LightSpeed'); 
fc = 300e6; 
mylambda = c/fc;
azangles = -90:0.5:90;
transmit_antennas = [16, 64];
Pt = 40;
L = length(transmit_antennas);

% --------- Plot 1: Single User MISO Beamforming --------------------------%

% --- Simulation Parameters -------%
theta = 0;
nUsers = 1;
fileLocation = '..Data/Beamforoming_investigation/single_user_mrt_beampattern.pdf';

% ---- Generate the MTR Precoders
w_single_mrt = generate_MRT_precoding_weights(nUsers, L, transmit_antennas, mylambda, theta, fc, Pt);
w_mrt = w_single_mrt;

% ----- Plot the resulting Beampattern 
plot_figure(L, transmit_antennas, w_mrt, azangles, fc, fileLocation, mylambda, theta);

% --------- Plot 2: Multi User MISO Beamforming --------------------------%

% --- Simulation Parameters -------%
theta = [-40, 0, 40];                       )
nUsers = length(theta);
fileLocation = 'Data/Beamforoming_investigation/multi_user_mrt_beampattern.pdf';

% ---- Generate the MTR Precoders
w_mrt = generate_MRT_precoding_weights(nUsers, L, transmit_antennas, mylambda, theta, fc, Pt);

% ----- Plot the resulting Beampattern 
plot_figure(L, transmit_antennas, w_mrt, azangles, fc, fileLocation, mylambda, theta);


% -------- Functions -------------------%

% ----- Generate the MRT Precoders ------------------%
function w_mrt = generate_MRT_precoding_weights(nUsers, L, transmit_antennas, mylambda, theta, fc, Pt)
    for idx = 1:L
        nTx = transmit_antennas(idx);
        
        % ---- MATLAB Phased Array Toolbox
        antennaArray = phased.ULA('NumElements', nTx, 'ElementSpacing', mylambda/2);
        delay = phased.ElementDelay('SensorArray', antennaArray);
        
        h = zeros(nTx, nUsers);
        for k = 1:nUsers
            tau = delay([theta(k); 0]);
            h(:, k) = exp(-1j*2*pi*fc*tau(:));
        end
        
        Px = Pt / nUsers;
        for k = 1:nUsers
            w_mrt(1:nTx, k, idx) = sqrt(Px) * h(:, k) / norm(h(:, k));
        end
    end
end

% ---------- Plot the Resulting Beampatterns -------------------%
function plot_figure(L, transmit_antennas, w_mrt, azangles, fc, fileLocation, mylambda, theta)
    figure
    AxesH = axes;
    hold on;
    labels = {};
    
    for idx = 1:L
        nTx = transmit_antennas(idx);
        
        antennaArray = phased.ULA('NumElements', nTx, 'ElementSpacing', mylambda/2);
        delay = phased.ElementDelay('SensorArray', antennaArray);
        
        W_k = w_mrt(1:nTx, :, idx);
        
        P_theta = zeros(1, length(azangles));
        for a = 1:length(azangles)
            tau = delay([azangles(a); 0]);
            sv  = exp(-1j*2*pi*fc*tau(:));
            P_theta(a) = sum(abs(W_k' * sv).^2);
        end
        
        P_dB = 10*log10(P_theta / max(P_theta));
        
        plot(azangles, P_dB, 'LineWidth', 1.4);
        labels{end+1} = sprintf('MRT, $M = %d$', nTx);
    end
    
    xline(theta, ':k', 'HandleVisibility', 'off');
    yline(-20, ':', '$-20$ dB', 'Interpreter', 'latex', 'FontSize', 12);
    
    xlabel('Azimuth angle $\theta$ (degrees)', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('Normalised radiated power (dB)', 'Interpreter', 'latex', 'FontSize', 12)
    legend(labels, 'Interpreter', 'latex', 'Location', 'south', 'FontSize', 10)
    
    xticks(-90:30:90);
    xaxisproperties = get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
    yaxisproperties = get(gca, 'YAxis');
    yaxisproperties.TickLabelInterpreter = 'latex';
    
    xlim([-90 90]);
    ylim([-40 2]);
    grid on;
    box on;
    
    InSet = get(AxesH, 'TightInset');
    set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3)-0.01, 1-InSet(2)-InSet(4)])
    exportgraphics(gcf, fileLocation, 'ContentType', 'vector', 'Resolution', 600)
end
