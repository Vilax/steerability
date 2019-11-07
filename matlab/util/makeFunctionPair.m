function [f1, f2, u1, u2, phi] = makeFunctionPair(alpha, numSamples)
    
    N = numel(alpha);
    u1 = zeros(size(alpha));
    u2 = zeros(size(alpha));
    u1(1:2:N) = alpha(1:2:N);
    u2(2:2:N) = alpha(2:2:N);
    u2 = [u2(2:end); 0] / sum(u2(:));
    u1 = [u1;0] / sum(u2(:));
    
    dphi = pi/(numSamples - 1);
    phi = [0:dphi:pi]';
    cosPhi = cos(phi);
    
    f1 = polyval(u1, cosPhi);
    f2 = polyval(u2, cosPhi);
end

