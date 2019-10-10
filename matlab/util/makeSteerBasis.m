function steerFiltVol = makeSteerBasis(filt, M)
% 
% makes steering basis given initial filter values and M number of basis
% filters to use. 
    angles = getSteerAngles2D(M);
    anglesDeg = angles * 180 / pi;
    [nrows, ncols] = size(filt);
    nfilts = M;
    steerFiltVol = zeros([nrows, ncols, nfilts]);
    
    for id = 1:nfilts
        angleRotate = anglesDeg(id);
        thisFilt = imrotate(filt, angleRotate, 'bilinear', 'crop');
        steerFiltVol(:,:,id) = thisFilt;
    end
end

