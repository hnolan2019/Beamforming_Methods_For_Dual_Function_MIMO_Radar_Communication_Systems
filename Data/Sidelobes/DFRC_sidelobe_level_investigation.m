% DFRC Sidelobe investigation and plotting
% Sidelobe Power vs SINR Threshold
clear; close all; clc;

% --- Assign colours following--- %
% Colours taken from https://rgbcolorpicker.com/
c_K2 = [0, 114,189]/255;
c_K4 = [217, 83,25]/255;
c_K6 = [237, 177, 32]/255;

% --------- Declare the three plot names --------%
fileLocation = 'DFRC_Sidelobe_Plot.pdf';

% --- Load combined sidelobe data ---%
data = readtable('updated_sidelobe_power_results.csv');

% --- Reshape into per-K columns for each method ---%
gamma = unique(data.Gamma_dB);
mse.K_2 = data.MSE_Sidelobe_Power_dB(data.K == 2);
mse.K_4 = data.MSE_Sidelobe_Power_dB(data.K == 4);
mse.K_6 = data.MSE_Sidelobe_Power_dB(data.K == 6);

maxmin.K_2 = data.MaxMin_Sidelobe_Power_dB(data.K == 2);
maxmin.K_4 = data.MaxMin_Sidelobe_Power_dB(data.K == 4);
maxmin.K_6 = data.MaxMin_Sidelobe_Power_dB(data.K == 6);

sca.K_2 = data.SCA_Sidelobe_Power_dB(data.K == 2);
sca.K_4 = data.SCA_Sidelobe_Power_dB(data.K == 4);
sca.K_6 = data.SCA_Sidelobe_Power_dB(data.K == 6);

% -------------Plot ------------------------%
figure;
AxesH = axes;
hold on; box on; grid on;

% --- Plot MSE --------- ---%
plot(gamma, mse.K_2, '-o', 'Color', c_K2,   'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin adaptive $\kappa$, $K=2$');
plot(gamma, mse.K_4, '-o', 'Color', c_K4, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin adaptive $\kappa$, $K=4$');
plot(gamma, mse.K_6, '-o', 'Color', c_K6, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin adaptive $\kappa$, $K=6$');

% --- MaxMin (fixed kappa) ---%
plot(gamma, maxmin.K_2, '--s', 'Color', c_K2,   'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin fixed $\kappa$, $K=2$');
plot(gamma, maxmin.K_4, '--s', 'Color', c_K4, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin fixed $\kappa$, $K=4$');

% --- Maxmin (adaptive kappa) ---%
plot(gamma, sca.K_2, ':^', 'Color', c_K2,   'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin adaptive $\kappa$, $K=2$');
plot(gamma, sca.K_4, ':^', 'Color', c_K4, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin adaptive $\kappa$, $K=4$');
plot(gamma, sca.K_6, ':^', 'Color', c_K6, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'MaxMin adaptive $\kappa$, $K=6$');

% Separate the legend - as it takes up a lot of space
% Method for two legends adapted from https://uk.mathworks.com/matlabcentral/answers/430791-how-to-add-a-second-legend-box-to-a-figure-without-new-plots

%----------- legend 1- The 3 Beamforming Methods -------- %
h_maxmin = plot(nan, nan, '--s',  'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4);
h_mse = plot(nan, nan, '-o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4);
h_sca = plot(nan, nan, ':^',  'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4);

% ----------- legend 2- The Number of Users -------- %
h_K2 = plot(nan, nan, '-', 'Color', c_K2, 'LineWidth', 3);
h_K4 = plot(nan, nan, '-', 'Color', c_K4, 'LineWidth', 3);
h_K6 = plot(nan, nan, '-', 'Color', c_K6, 'LineWidth', 3);

hold off;

% --- Axes formatting --- ----%
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 0.8;
ax.TickDir = 'in';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;

ax.XLim = [gamma(1), gamma(end)];
ax.YLim = [2, 9];
ax.XTick = gamma(1):2:gamma(end);          
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

xlabel('SINR threshold $\Gamma$ (dB)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Average sidelobe power $P_{side}$ (dB)',  'Interpreter', 'latex', 'FontSize', 12);

% --- Place main axes tight ---%
%
% Method adapted from https://uk.mathworks.com/matlabcentral/answers/471605-how-to-use-position-property-to-make-axes-height-tight-to-figure-height
InSet = ax.TightInset;
ax.Position = [InSet(1:2), 1-InSet(1)-InSet(3)-0.01, 1-InSet(2)-InSet(4)];

% --- Legend 1: Method ---%
% The code is adpated from https://blogs.mathworks.com/pick/2011/06/24/flexible-legends/
leg1 = legend(ax, [h_mse, h_maxmin, h_sca], ...
              {'MSE','MaxMin (fixed $\kappa$)', 'MaxMin (adaptive $\kappa$)'}, ...
              'Interpreter', 'latex', 'FontSize', 9, ...
              'Box', 'on', 'EdgeColor', [0.5 0.5 0.5], ...
              'AutoUpdate', 'off');

leg1.Position(1) = ax.Position(1) + 0.02;
leg1.Position(2) = ax.Position(2) + ax.Position(4) - leg1.Position(4) - 0.02;

% --- Legend 2: Number of users ---%
ax2 = axes('Position', ax.Position, 'Visible', 'off');
leg2 = legend(ax2, [h_K2, h_K4, h_K6], ...
              {'$K = 2$', '$K = 4$', '$K = 6$'}, ...
              'Interpreter', 'latex', 'FontSize', 9, ...
              'Box', 'on', 'EdgeColor', [0.5 0.5 0.5], ...
              'AutoUpdate', 'off');

% Place directly below legend 1
leg2.Position(1) = leg1.Position(1);
leg2.Position(2) = leg1.Position(2) - leg2.Position(4);

hold off;

% --- Export to PDF ---%
exportgraphics(gcf, fileLocation, 'ContentType', 'vector', 'Resolution', 600);