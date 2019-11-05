function dirAngles = getDirectionAngles(N, axis)

    assert(norm(axis) ~= 0)
    
    alpha = axis(1) /norm(axis);
    beta = axis(2)/norm(axis);
    gamma = axis(3)/norm(axis);
    
    powers = dirCosPowers3D(N);
    alphaPowers = powers(:,1);
    betaPowers = powers(:,2);
    gammaPowers = powers(:,3);
    
    dirAlpha = (alpha .^ alphaPowers);
    dirBeta = (beta .^ betaPowers);
    dirGamma = (gamma .^ gammaPowers);
    dirAngles = dirAlpha .* dirBeta .* dirGamma;    
end

