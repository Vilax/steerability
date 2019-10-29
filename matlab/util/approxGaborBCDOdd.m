function [x, approximator, gaborSigma] = ...
                    approxGaborBCDOdd(n, sigmaDenom, y0Denom, nloops)
% APPROXGABORBCD    approximate odd Gabor Fourier Transform using BCD
% creates a initial gabor based on input, and then uses BCD on - in order,
% coefficients of polynomial of the radius, coefficients of polynomial of
% the angle, sigma size of gabor approximator, sigma size of gabor to be
% approximated
    
    DEFAULT_NLOOPS = 10; 
    if nargin < 4
            nloops = DEFAULT_NLOOPS;
    end
    
    % first initialize gabor
    gaborSigma = n/sigmaDenom;
    y0 = n/y0Denom;
    [X, Y] = meshgrid(-n:n);
    R = sqrt(X.^2 + Y.^2); 
    
    % build gabor fourier transform
    gaborFT = 0.5 * exp( -((X).^2+(Y-y0).^2)/gaborSigma^2) - ...
        0.5 * exp( -((X).^2+(Y+y0).^2)/gaborSigma^2);
    
    angle = atan2(X,Y);
    Theta = angle .* (angle >= 0) - (angle .* (angle < 0));
    cosTheta = cos(Theta);
    
    % now create initial set of values for the approximator
    initFun = @(x) sum(sum((gaborFT - (exp(-((R-y0)/(n/x(1))).^2)) .* ...
                (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
                (polyval([x(6), 0, x(5), 0], cosTheta))).^2, 2), 1);
    xinit = [sigmaDenom, 1, 0, 0, 0, 0];
    [x, fval] = fminsearch(initFun, xinit); 
    
    for iloop = 1:nloops
        sigma = x(1);
        b = [x(5), x(6)];
        polyRfunc = @(a) sum(sum((gaborFT-(exp(-((R-y0)/(n/sigma)).^2)) ...
                 .* (polyval([a(3), a(2), a(1)], (R-y0))) .* ...
                (polyval([b(2), 0, b(1), 0], cosTheta))).^2, 2), 1);
        ainit = [x(1), x(2), x(3)];
        [a, fval] = fminsearch(polyRfunc, ainit);
        % collect optimum values
        
        x = [sigma, a, b];
        assert(numel(x) == 6)
        
        polyAnglefunc = @(b) sum(sum((gaborFT-(exp(-((R-y0)/(n/sigma)).^2))...
                .* (polyval([a(3), a(2), a(1)], (R-y0))) .* ...
                (polyval([b(2), 0, b(1), 0], cosTheta))).^2, 2), 1);
        binit = [x(5), x(6)];
        [b, fval] = fminsearch(polyAnglefunc, binit);

        x = [sigma, a, b];
        assert(numel(x) == 6)
        approxSigmafunc = @(sigma) sum(sum((gaborFT - ...
                    (exp(-((R-y0)/(n/sigma)))) .* ...
                    (polyval([a(3), a(2), a(1)], (R-y0))) .* ...
                    (polyval([b(2), 0, b(1), 0], cosTheta))).^2, 2), 1);
        sigmaInit = x(1);
        [sigma, fval] = fminsearch(approxSigmafunc, sigmaInit);

        x = [sigma, a, b];
        assert(numel(x) == 6)
        approximator = (exp(-((R-y0)/(n/sigma)).^2)) .* ...
                    (polyval([a(3), a(2), a(1)], (R-y0))) .* ...
                    (polyval([b(2), 0, b(1), 0], cosTheta));  
        gaborAdjustSigmafunc = @(x) sum(sum( (0.5 * ...
                (exp( -((X).^2+(Y-y0).^2)/(x(1))^2) - ...
                exp(-((X).^2+(Y+y0).^2)/(x(1))^2)) - approximator).^2, 1)...
                    , 2);  
        xinit = gaborSigma;
        [x, fval] = fminsearch(gaborAdjustSigmafunc, xinit);
        
        gaborSigma = x; 
        x = [sigma, a(1), a(2), a(3), b(1), b(2)];
        assert(numel(x) == 6)
        fval
    end
end

