% build gabor

n = 80;
[X, Y, Z] = meshgrid(-n:n);
R = sqrt(X.^2 + Y.^2 + Z.^2);

sigma = n/6;
z0 = n/4;
gaborEvenFT3D = 0.5 * (exp( -(X.^2 + Y.^2 + (Z-z0).^2)/sigma^2)) + ...
                (exp( -(X.^2 + Y.^2 + (Z+z0).^2)/sigma^2));

gaborOddFT3D = 0.5 * (exp( -(X.^2 + Y.^2 + (Z+z0).^2)/sigma^2)) - ...
                (exp( -(X.^2 + Y.^2 + (Z-z0).^2)/sigma^2));      
            
[u,v,w] = project2Sphere(X,Y,Z);

cosPhi = w;

fun = @(x)sum(sum(sum((gaborEvenFT3D - (exp(-(R/(n/x(1))).^2)) .* ... 
            (polyval([x(4),x(3),x(2)],R)) .* ...
            (polyval([x(6),0, x(5),0,0], cosPhi)) ).^2 ,1),2),3)
        
xinit = [5, 2, 2, 0.7, 1, 1.3];        
[x, fval] = fminsearch(fun, xinit);
current_x = x;
current_fval = fval;

x1range = [3:8];
polyxvalrange = [0.2:0.4:3];

for idx1 = 1:numel(x1range)
    x1 = polyxvalrange(idx1);
    for idx2 = 1:numel(polyxvalrange)
        x2 = polyxvalrange(idx2);
        for idx3 = 1:numel(polyxvalrange)
            x3 = polyxvalrange(idx3);
            for idx4 = 1:numel(polyxvalrange)
                x4 = polyxvalrange(idx4);
                for idx5 = 1:numel(polyxvalrange)
                    x5 = polyxvalrange(idx5);
                    for idx6 = 1:numel(polyxvalrange)
                        x6 = polyxvalrange(idx6);
                        
                        xinit = [x1, x2, x3, x4, x5, x6];
                        [x, fval] = fminsearch(fun,xinit);
                        if fval < current_fval
                            current_fval = fval; 
                            current_x = x;
                        end
                    end
                end
            end
        end
    end
end
x = current_x
approximator = (exp(-(R/(n/x(1))).^2)) .* ... 
            (polyval([x(4),x(3),x(2)],R)) .* ...
            (polyval([x(6),0, x(5),0,0], cosPhi));
        
WriteMRC(gaborEvenFT3D, 1, 'gaborFiltEven3D.mrc')
WriteMRC(approximator, 1, 'gaborApproximatorEven3D.mrc')

% here begins attempts to approximate with a shift term in approximator
% gaussian term
n = 40;
[X, Y, Z] = meshgrid(-n:n);
R = sqrt(X.^2 + Y.^2 + Z.^2);

sigma = n/6;
z0 = n/4;
gaborEvenFT3D = 0.5 * (exp( -(X.^2 + Y.^2 + (Z-z0).^2)/sigma^2)) + ...
                (exp( -(X.^2 + Y.^2 + (Z+z0).^2)/sigma^2));

gaborOddFT3D = 0.5 * (exp( -(X.^2 + Y.^2 + (Z+z0).^2)/sigma^2)) - ...
                (exp( -(X.^2 + Y.^2 + (Z-z0).^2)/sigma^2));      
            
[u,v,w] = project2Sphere(X,Y,Z);

cosPhi = w;

func = @(x) sum(sum(sum((gaborEvenFT3D - (exp(-((R-z0)/(n/x(1))).^2)) .* ... 
            (polyval([x(4),x(3),x(2)],(R-z0))) .* ...
            (polyval([x(6),0, x(5),0,0], cosPhi))).^2, 1), 2), 3)
        
xinit = [sigma, 0, 0 ,1,0,0];
options = optimset('PlotFcns',@optimplotfval);
[x,fval] = fminsearch(func, xinit, options);


    approximator = (exp(-((R-z0)/(n/x(1))).^2)) .* ... 
                (polyval([x(4),x(3),x(2)],(R-z0))) .* ...
                (polyval([x(6),0, x(5),0,0], cosPhi));

    WriteMRC(gaborEvenFT3D, 1, 'gaborFiltEven3D.mrc')
    WriteMRC(approximator, 1, 'gaborApproximatorEven3D.mrc')