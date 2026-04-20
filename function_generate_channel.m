% --------------- The Wireless Channel Model --------------- %
function h = function_generate_channel(M, K)


    % ---- Function to generate the communication channel -------------%
    %
    % Inputs - number of transmit antennas and number of users
    % Outputs - the channel vector h

    h = zeros(M, K);                        % Allocate memory for channel
    large_fading_scale = zeros(1, K);       % Allocate memory for (LSF) array    
    min_dist = 0.05;                        % Minimum distance in km
    max_dist = 2;                           % Maximum distance in km
    shadowdev = 10;                         % Shadowning value in dB
    
    % Large Scale Fading value for each user at distance d from Tx
    for iUsers = 1:length(large_fading_scale)
        user_k_d = min_dist+rand()*(max_dist-min_dist);
        pathloss_k = -(128.1+37.6*log10(user_k_d));
        betadB = pathloss_k + shadowdev*randn(1);
        large_fading_scale(iUsers) = 10^(betadB/10);
    end
    
    % Large and Small Scale fading factors
    alpha = (randn(M, K)+1i*randn(M, K))/sqrt(2);
    beta = sqrt(large_fading_scale);                
    
    % ------- Channel Noise -----------------------%
    NF = 9;                                 % Noise factor (dB)
    No = -174 - 30 ;                        % dBW/Hz
    B = 10*log10(10e6);                     % 10 MHz in dB
    noisepowerdB = NF+No+B;
    noisepower=10^(noisepowerdB/10);        % Noise power in linear scale
    
    % The channel with fading and adjusted for noise
    h = alpha*diag(beta)/sqrt(noisepower);
end
