function dirVals = computeSteerDirection2D(N, Theta, nonzeroCoeff)
    
    if nargin < 3
        nonzeroCoeff = ones([2*N+1,1]);
    end
    numNonzero = validateNonzeroCoefficients(nonzeroCoeff, N);
    powersBase = [-N:N]';
    powerSpec = powersBase(nonzeroCoeff == 1);
    
    dirVals = exp(1i*(Theta .* powerSpec));

end

