function normVals = computeNormalizationConstant2D(N)
% COMPUTENORMALIZATIONCONSTANT2D    Compute normalization constants for 2d
% s1 basis
% normVals=COMPUTENORMALIZATIONCONSTANT2D(N) computes for up
% to bandwidth N the normalization constants for the basis vectors 1, x^2,
% x^4,... x^N by saving the integral from 0 to pi/2 of cos^{2m} xdx where m
% ranges from 0 to N, where N=2M is assumed even.
% result will be an array normVals where normVals(k) = 1/(int_0^pi/2
% cos^{4(k-1)}x dx), that is, writing alpha_k = normVals(k+1), that is, 
% adjusting for 1-indexing in matlab, alpha_k is the
% constant s.t. each alpha_k is the coefficient to x^{2k}, for k=0,..,M,
% such that the int_0^pi/2 (alpha_k x^2k)^2 dx is 1.

    assert(N>=1);
    N = validateBandwidth(N);
    M = N/2;
    integralVals = zeros([M+1, 1]);
    
    % initialize first steps
    integralVals(1) = pi/2;
    
    % go through array
    % note that in general int_0^pi/2 cos^m xdx = (m-1)/m * int_0^pi/2
    % cos^{m-2} x dx = (m-1)(m-3)/m(m-2) int_0^pi/2 cos^{m-4}xdx
    for ival = 2:numel(integralVals)
        m = 4*(ival-1); % retrieve power
        scaling = (m-1)*(m-3)/(m*(m-2));
        integralVals(ival) = scaling* integralVals(ival-1);    
    end
    normVals = sqrt(integralVals);
end

