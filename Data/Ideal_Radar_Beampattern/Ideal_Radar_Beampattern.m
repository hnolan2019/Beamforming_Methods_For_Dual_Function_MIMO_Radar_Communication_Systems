% Adapted from MATLAB - help cener 
% https://uk.mathworks.com/help/phased/ug/waveform-design-for-a-dual-function-mimo-radcom-system.html

tgtAz = [-60 0 60];                 % Azimuths of the targets of interest
tgtRng = [5.31e3, 6.23e3, 5.7e3];   % Ranges of the targets of interest

ang = linspace(-90, 90, 200);       % Grid of azimuth angles
beamwidth = 20;                     % Desired beamwidth

% Desired beam pattern
idx = false(size(ang));
for i = 1:numel(tgtAz)
    idx = idx | ang >= tgtAz(i)-beamwidth/2 & ang <= tgtAz(i)+beamwidth/2;
end

Bdes = zeros(size(ang));
Bdes(idx) = 1;
blue_c  = [21,0,245]/255; 

% --- Plot --- %
figure;
AxesH = axes;
hold on; box on; grid on;

fileLocation = 'ideal_radar_beampattern.pdf';

plot(ang, Bdes , 'Color', blue_c,   'LineWidth', 2, 'HandleVisibility', 'off');
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
ax.YLim = [0, 1.1];
ax.XTick = -90:30:90;

ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';

xlabel('Azimuth angle $\theta$ (deg)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Normalised Spatial Response',  'Interpreter', 'latex' ,'FontSize', 12)

% --- Place main axis --- %
InSet = ax.TightInset;
ax.Position = [InSet(1:2), 1-InSet(1)-InSet(3)-0.01, 1-InSet(2)-InSet(4)];

% --- Export to PDF --- %
exportgraphics(gcf, fileLocation, 'ContentType', 'vector', 'Resolution', 600);