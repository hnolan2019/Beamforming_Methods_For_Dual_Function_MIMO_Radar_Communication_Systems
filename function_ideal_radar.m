% -----------------Ideal Radar beam pattern ------------------%
function des = function_ideal_radar(ang, tgtAz, beamwidth)

% ---- Function to plot the ideal radar waveform -------------%
%
% Inputs - grid points, target angles and the beamwdith
% Outputs - the ideal radar beampattern

    idx = false(size(ang));
    for i = 1:length(tgtAz)
        idx = idx | ang >= tgtAz(i)-beamwidth/2 & ang <= tgtAz(i)+beamwidth/2;
    end
    des = zeros(size(ang));
    des(idx) = 1;
end