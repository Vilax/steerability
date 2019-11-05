function polMat = cartMat2PolMat(cartMat, dims)
% Given a matrix cartMat, converts to a matrix polMat whose column
% coordinate is theta and whose row coordiante is r. Theta here is angle to
% north pole y-axis

    xmin = dims(1);
    xmax = dims(2);
    ymin = dims(3);
    ymax = dims(4);
    
    [X,Y] = meshgrid([xmin:xmax], [ymin:ymax]);
    
    R = sqrt(X.^2 + Y.^2);
    angle = atan2(X,Y);
    Theta = angle .* (angle >= 0) - (angle .* (angle < 0));
    Theta = Theta * 180 / pi;
    
    Rqmax = floor(max(R(:)));
    Rqrange = [0:Rqmax];
    
    thetaRange = [0:180];
    
    [Thetaq, Rq] = meshgrid(thetaRange, Rqrange);
    
    Vq = interp2(Theta(:), R(:), cartMat(:), Thetaq, Rq);
    
    poltMat = Vq;
    
    
end

