function filt = makeSteerFilt3D(n, r0, sigmaDenom, constructMethod, varargin)
% MAKESTEERFILT3D
% 
% filt=MAKESTEERFILT3D(n,r0,sigmaDenom,'poly', p) 
% filt=MAKESTEERFILT3D(n,r0,sigmaDenom,'interpolation',f,phi) expects two
% arguments after the construct method which gives the function values and
% corresponding phi (angle to positive z-axis) values (expected to range
% from 0 to pi)

    [x,y,z] = meshgrid(-n:n, -n:n, -n:n);
    r = sqrt(x.^2+y.^2+z.^2);
    sigma = n/sigmaDenom;
    
    g = exp(-(r-r0).^2/(2*sigma^2));
    
    [u, v, w] = project2Sphere(x,y,z);
     df
    if isequal(constructMethod, 'poly')
        assert(numel(varargin) == 1);
        p = varargin{1};
        sphericalVals = polyval(p,w);
    elseif isequal(constructMethod, 'interpolation')
        assert(numel(varargin) == 2);
        f = varargin{1};
        phi = varargin{2};
        angles = acos(w);
        a = reshape(angles, [numel(angles(:)) 1]);
        values = interp1(phi, f, a); % should this use a? test this
        sphericalVals = reshape(values, size(angles));
    end
    
    filt = g .* sphericalVals;
    
end

