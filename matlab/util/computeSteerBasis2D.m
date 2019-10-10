function steerBasisVals = computeSteerBasis2D(N, nonzeroCoeff)
% assumes nonzeroCoeff is a vector of bools
    numNonzero = validateNonzeroCoefficients(nonzeroCoeff, N);
    powersBase = [-N:N]';
    powerSpec = powersBase(nonzeroCoeff == 1);
    steerTheta = getSteerAngles2D(numNonzero);
    steerThetaBase = repmat(steerTheta, [numNonzero,1]);
    steerBasisVals = exp(1i*(steerThetaBase .* powerSpec));
end

