% Transmit Beampattern Plot
clear; close all; clc;

% --------- Declare the four plot names --------%
fileLocation1 = 'transmit_beampattern_2_users_10dB.pdf';
fileLocation2 = 'transmit_beampattern_4_users_10dB.pdf';
fileLocation3 = 'transmit_beampattern_6_users_10dB.pdf';
fileLocation4 = 'transmit_beampattern_2_users_2dB.pdf';

% --------- Read in the 3 data tables ----------%
D1= readtable('updated_beampattern_results_k_2aB_10.csv');  % 2 users 10dB
D2 = readtable('beampattern_results_k_4_db_10.csv');        % 4 users 10dB
D3 = readtable('beampattern_results_k_6_dB_10.csv');        % 6 users 10dB
D4 = readtable('beampattern_results_2_users_2dB.csv');      % 6 users 10dB

% ----- Plot 1: 2 Users 10dB ------%
az1 = D1.Azimuth_deg;
plot_transmit_beampattern(D1, az1, fileLocation1);

% ------ Plot 2: 4 Users 10dB ----------------------%
az2 = D2.Azimuth_deg;
plot_transmit_beampattern(D2, az2, fileLocation2);

% ------ Plot 3: 6 Users 10dB ----------------------%
az3 = D3.Azimuth_deg;
plot_transmit_beampattern(D3, az3, fileLocation3);

% ------ Plot 4: 2 Users 2dB ----------------------%
az4 = D4.Azimuth_deg;
plot_transmit_beampattern(D4, az4, fileLocation4);

% ------- Convert to db-------------------------%
function dB = convert_to_db(x)
    dB = 10 * log10(max(x, 1e-12));
end

function data = check_for_NaNs(datefile)

    % MaxMin - most likely to fail
    valid = ~isnan(datefile.MaxMin);
    
    data.MSE = datefile.MSE(valid);
    data.FinalMaxMin = datefile.FinalMaxMin(valid);
    data.MaxMin = datefile.MaxMin(valid);
    data.RadarOnly = datefile.RadarOnly(valid);
end

function plot_transmit_beampattern(D, az, fileLocation)

    % --- Assign colours following--- %
    % Colours taken from https://rgbcolorpicker.com/
    c_radar = [0,0,0]/255;                      % black 
    c_mse = [242,10, 10]/255;                   % blue
    c_maxmin = [21,0,245]/255;                  % red
    c_finalmm = [21,0,245]/255;                 % red

    % --- Plot --- %
    figure;
    AxesH = axes;
    hold on; box on; grid on;
    
    plot(az, convert_to_db(D.RadarOnly),  '-',  'Color', c_radar,   'LineWidth', 1.4, 'HandleVisibility', 'off');
    plot(az, convert_to_db(D.MSE), '--', 'Color', c_mse, 'LineWidth', 1.4, 'HandleVisibility', 'off');
    plot(az, convert_to_db(D.MaxMin), ':',  'Color', c_maxmin,  'LineWidth', 1.6, 'HandleVisibility', 'off');
    plot(az, convert_to_db(D.FinalMaxMin), '-.', 'Color', c_finalmm, 'LineWidth', 1.4, 'HandleVisibility', 'off');
    
    % --- legen --- %
    h_radar = plot(nan, nan, '-',  'Color', c_radar,'LineWidth', 1.4);
    h_mse = plot(nan, nan, '--', 'Color', c_mse, 'LineWidth', 1.4);
    h_maxmin = plot(nan, nan, ':',  'Color', c_maxmin, 'LineWidth', 1.6);
    h_finalmm = plot(nan, nan, '-.', 'Color', c_finalmm, 'LineWidth', 1.4);
    
    hold off;
    
    % --- Axes formatting --- %
    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 0.8;
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.5;
    
    ax.XLim  = [-90, 90];
    ax.YLim = [-10, 20];
    ax.XTick = -90:30:90;
    
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    
    xlabel('Azimuth angle $\theta$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Transmit beampattern $P(\theta; \mathbf{R})$ (dB)', 'Interpreter', 'latex', 'FontSize', 12);
    
    % --- Place main axis --- %
    % Method for two legends adapted from https://uk.mathworks.com/matlabcentral/answers/430791-how-to-add-a-second-legend-box-to-a-figure-without-new-plots
    InSet = ax.TightInset;
    ax.Position = [InSet(1:2), 1-InSet(1)-InSet(3)-0.01, 1-InSet(2)-InSet(4)];
    
    % --- Legend --- %
    % The code is adpated from https://blogs.mathworks.com/pick/2011/06/24/flexible-legends/
    leg = legend(ax, [h_radar, h_mse, h_maxmin, h_finalmm], ...
        {'Radar only', 'MSE', 'MaxMin (fixed $\kappa$)', 'MaxMin (adaptive $\kappa$)'}, ...
        'Interpreter', 'latex', 'FontSize', 9, ...
        'Box', 'on', 'EdgeColor', [0.5 0.5 0.5], ...
        'AutoUpdate', 'off');
    leg.Position(1) = ax.Position(1) + ax.Position(3) - leg.Position(3) - 0.01;
    leg.Position(2) = ax.Position(2) + 0.01; 
    
    % --- Export to PDF --- %
    exportgraphics(gcf, fileLocation, 'ContentType', 'vector', 'Resolution', 600);
end

