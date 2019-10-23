function filt = makeSteerFilt(n, f, theta,r0,sigmaDenom)
% MAKESTEERFILT 
%
    [x,y]=meshgrid(-n:n,-n:n);
    r=sqrt(x.^2+y.^2);
    sigma=n/sigmaDenom;

    g=exp(-(r-r0).^2/(2*sigma^2));
    %g = ones(size(g))
    angle=atan2(x,y);
    % angle to north pole must be in range 0 to pi
    angle = angle .* (angle >= 0) - (angle .* (angle < 0));
    angleSize=size(angle);
    a=reshape(angle,[prod(angleSize) 1]);
    values=interp1(theta,f,a);
    
    filt=g.*reshape(values,angleSize);
end

% function filt = makeSteerFilt(n, f, theta,r0)
% % MAKESTEERFILT 
% %
%     [x,y]=meshgrid(-n:n,-n:n);
%     r=sqrt(x.^2+y.^2);
%     sigma=n/6;
% 
%     g=exp(-(r-r0).^2/(2*sigma^2));
%     angle=abs(atan2(x,y));
%     angle=angle.*(angle<=pi/2)+(pi-angle).*(angle>pi/2);
%     angleSize=size(angle);
%     a=reshape(angle,[prod(angleSize) 1]);
%     values=interp1(theta,f,a);
%     
%     filt=g.*reshape(values,angleSize);
%     
% end