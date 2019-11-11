function filt = makeRotatedFilt(n,f,theta,r0,sigmaDenom,angleRotate)

    [x,y]=meshgrid(-n:n,-n:n);

    %g = ones(size(g))
    angle=atan2(x,y);
    % range now from -pi to pi
    angleRotate = mod(angleRotate, 360);
    angle = angle - (pi/180)*angleRotate;
    angle = mod(angle, 2*pi);
    
    % angle to north pole must be in range 0 to pi
%     angle = angle .* (angle >= 0) - (angle .* (angle < 0));
    angle = angle .* (angle <= pi) + ((2*pi - angle) .* (angle > pi));
    angleSize=size(angle);
    a=reshape(angle,[prod(angleSize) 1]);
    values=interp1(theta,f,a);
    
    g = makeRadialFunction('gaussian', n,r0,sigmaDenom);
    g = makeRadialFunction('spline', n, r0);
    filt=g.*reshape(values,angleSize);
end

