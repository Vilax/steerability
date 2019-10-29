function [x, approximator, gaborSigma, gaborFT] = ...
                    approxGaborBCD(n, sigmaDenom, y0Denom, nloops, parity)
% APPROXGABORBCD    approximate Gabor Fourier Transform using BCD
% creates a initial gabor based on input, and then uses BCD on - in order,
% coefficients of polynomial of the radius, coefficients of polynomial of
% the angle, sigma size of gabor approximator, sigma size of gabor to be
% approximated
    DEFAULT_PARITY = 'odd';
    DEFAULT_NLOOPS = 10; 
    if nargin < 5
        parity = DEFAULT_PARITY;
        if nargin < 4
            nloops = DEFAULT_NLOOPS;
        end
    end
    
    % first initialize gabor
    if isequal(parity, 'even')
        [x, approximator, gaborSigma] = approxGaborBCDEven(n,sigmaDenom,...
                                    y0Denom, nloops);
        gaborFT = 0.5 * exp( -((X).^2+(Y-y0).^2)/gaborSigma^2) + ...
                    0.5 * exp( -((X).^2+(Y+y0).^2)/gaborSigma^2);
    else
        [x, approximator, gaborSigma] = approxGaborBCDOdd(n,sigmaDenom,...
                                    y0Denom, nloops);
        gaborFT = 0.5 * exp( -((X).^2+(Y-y0).^2)/gaborSigma^2) - ...
                    0.5 * exp( -((X).^2+(Y+y0).^2)/gaborSigma^2);
    end 
end

