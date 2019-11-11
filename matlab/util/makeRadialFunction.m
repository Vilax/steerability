function g = makeRadialFunction(type, varargin)
    if isequal(type, 'gaussian')
        n = varargin{1};
        r0 = varargin{2};
        sigmaDenom = varargin{3};
        g = makeGaussianRadial(n, r0, sigmaDenom);
    elseif isequal(type, 'spline')
        n = varargin{1};
        u = varargin{2};
        g = makeSplineRadial(n,u);
    end
end

function g=makeGaussianRadial(n,r0,sigmaDenom)
    [x,y] = meshgrid(-n:n, -n:n);
    r = sqrt(x.^2 + y.^2);
    sigma = n/sigmaDenom;
    g=exp(-(r-r0).^2/(2*sigma^2));
end

function g = makeSplineRadial(n, r0)
    [x,y] = meshgrid(-n:n, -n:n);
    r = sqrt(x.^2 + y.^2);   
    [theta,r]=cart2pol(x,y);
    
    r_n=1.5*r/r0;
    r_n=(r_n-1.5)+1.5;
    g=b_spline(r_n);
end

