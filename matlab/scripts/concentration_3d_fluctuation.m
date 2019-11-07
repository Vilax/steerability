capSizeDenom = [2.1:0.1:5]; 
minVals = zeros(size(capSizeDenom));
for iDenom = [1:numel(capSizeDenom)]
    capSize = pi/capSizeDenom(iDenom);
    [f,u,~,phi] = steer3dGeneral(capSize,3);
    midpoint = floor(numel(f)/2);
    tmp = f;
    tmp(midpoint:end) = abs(tmp(midpoint:end));
    minVals(iDenom) = min(tmp(:));
end

figure; plot(pi./capSizeDenom, minVals);