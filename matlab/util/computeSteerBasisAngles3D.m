function steerAngles = computeSteerBasisAngles3D(N, DirCos)

    M = (N+1)*(N+2)/2;

    powers = dirCosPowers3D(N);
    alphaPowers = powers(:,1);
    betaPowers = powers(:,2);
    gammaPowers = powers(:,3);
    
    steerAlphaBase = repmat(DirCos(:,1)', [M,1]);
    steerBetaBase = repmat(DirCos(:,2)', [M,1]);
    steerGammaBase = repmat(DirCos(:,3)', [M,1]);
    
    steerAlpha = (steerAlphaBase .^ alphaPowers);
    steerBeta = (steerBetaBase .^ betaPowers);
    steerGamma = (steerGammaBase .^ gammaPowers);
    steerAngles = steerAlpha .* steerBeta .* steerGamma;
end

