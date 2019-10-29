function approximator = buildGaborApproximator(x, n, y0, parity)
% BUILDGABORAPPROXIMATOR 

    DEFAULT_PARITY = 'odd';
    if nargin < 4
        parity = DEFAULT_PARITY;
    end
    
    [X, Y] = meshgrid(-n:n);
    R = sqrt(X.^2+Y.^2);
    angle = atan2(X,Y);
    Theta = angle .* (angle >= 0) - (angle .* (angle < 0));
    cosTheta = cos(Theta);
    
    if isequal(parity, 'odd')
        approximator = (exp(-((R-y0)/(n/x(1))).^2)) .* ...
                        (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
                        (polyval([x(6), 0, x(5), 0], cosTheta));
    else
        approximator = (exp(-((R-y0)/(n/x(1))).^2)) .* ...
                        (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
                        (polyval([x(6), 0, x(5), 0, 0], cosTheta));
    end
end

