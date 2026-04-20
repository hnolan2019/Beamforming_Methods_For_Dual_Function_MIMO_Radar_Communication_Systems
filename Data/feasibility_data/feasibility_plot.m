% Feasibility Comparison Plot
clear; close all; clc;

% --- Read in the data -------%
maxmin = readtable('feasibility_MaxMin.csv');
mse = readtable('feasibility_MSE.csv');
sca = readtable('feasibility_Proposed_SCA.csv');

%------ Output path ---%
fileLocation = 'feasibility_comparison.pdf';

gamma = maxmin.Gamma_dB;

% --- Assign colours following--- %
% Colours taken from https://rgbcolorpicker.com/
c_K2 = [0, 114,189]/255;
c_K4 = [217, 83,25]/255;
c_K6 = [237, 177, 32]/255;


% ------- Plot the Feasibility---
figure;
AxesH = axes;
hold on; box on; grid on;

% Plot Max-Min with fixed kappa
plot(gamma, maxmin.K_2/100, '--s', 'Color', c_K2, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');
plot(gamma, maxmin.K_4/100, '--s', 'Color', c_K4, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');
plot(gamma, maxmin.K_6/100, '--s', 'Color', c_K6, 'MarkerFaceColor', 'w',  'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');

% Plot MSE
plot(gamma, mse.K_2/100, '-o', 'Color', c_K2, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');
plot(gamma, mse.K_4/100, '-o', 'Color', c_K4, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');
plot(gamma, mse.K_6/100, '-o', 'Color', c_K6, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');

% Plot Max-Min with adjust kappa
plot(gamma, sca.K_2/100, ':^', 'Color', c_K2, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');
plot(gamma, sca.K_4/100, ':^', 'Color', c_K4, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');
plot(gamma, sca.K_6/100, ':^', 'Color', c_K6, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.4, 'HandleVisibility', 'off');


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
ax.YLim = [-0.02, 1.02];
ax.XTick = gamma(1):2:gamma(end);
ax.YTick = 0:0.2:1;

ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

xlabel('SINR threshold $\Gamma$ (dB)', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Feasibility probability', 'Interpreter', 'latex', 'FontSize', 13);

% ----- place main axis ------- %
InSet = ax.TightInset;
ax.Position = [InSet(1:2), 1-InSet(1)-InSet(3)-0.01, 1-InSet(2)-InSet(4)];

%-----legend 1 box
leg1 = legend(ax, [h_mse, h_maxmin, h_sca], { 'MSE', 'MaxMin (fixed $\kappa$)', 'MaxMin (adaptive $\kappa$)'}, 'Interpreter', 'latex', 'FontSize', 9, ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5], ...
    'AutoUpdate', 'off');
leg1.Position(1) = ax.Position(1) + ax.Position(3) - leg1.Position(3) - 0.01;
leg1.Position(2) = ax.Position(2) + ax.Position(4) - leg1.Position(4) - 0.01;

% --- Legend 2, number of users per plot
ax2 = axes('Position', ax.Position, 'Visible', 'off');
leg2 = legend(ax2, [h_K2, h_K4, h_K6], {'$K = 2$', '$K = 4$', '$K = 6$'}, ...
    'Interpreter', 'latex', 'FontSize', 9, ...
    'Box', 'on', 'EdgeColor', [0.5 0.5 0.5], ...
    'AutoUpdate', 'off');

% Put directly below 1
leg2.Position(1) = leg1.Position(1) + leg1.Position(3)-leg2.Position(3);
leg2.Position(2) = leg1.Position(2) - leg2.Position(4) - 0.005;

% export to pdf
exportgraphics(gcf, fileLocation, 'ContentType', 'vector', 'Resolution', 600);
