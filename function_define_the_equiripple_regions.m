% ------- Function to partition the grid into 3 different regions ------%

function [passband_idx, transition_idx, stopband_idx] = function_define_the_equiripple_regions(ang, beamwidth, delta_trans, tgtAz)

    % Inputs - grid points, beamwidth (passband), delta_trans(transisition band),
    %           tgtAz (target locations)
    %          
    % Outputs - grind index for passband, transition band and stopband

    passband_idx = false(size(ang));
    transition_idx = false(size(ang));
    
    for idx = 1:length(tgtAz)
        passband_idx = passband_idx | (ang >= tgtAz(idx) - beamwidth/2 & ang <= tgtAz(idx) + beamwidth/2);
        transition_idx = transition_idx | (ang >= tgtAz(idx) - beamwidth/2  - delta_trans & ang <= tgtAz(idx) + beamwidth/2  + delta_trans);
    end
    stopband_idx = ~transition_idx;
end