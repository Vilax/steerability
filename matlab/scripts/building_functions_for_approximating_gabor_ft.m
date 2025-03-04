% script building and demonstrating functions that use steerable filters to
% approximate the fourier transform of gabor functions
close all; clear all;
% first create a fourier transform of a gabor filter
% n = 100; % size of image
% [X, Y] = meshgrid(-n:n);
% R = sqrt(X.^2 + Y.^2);
% 
% sigma = n/6; 
% y0 = n/3;
% gaborEvenFT = 0.5 * (exp( -((X).^2+(Y-y0).^2)/sigma^2) + ...
%                     exp(-((X).^2+(Y+y0).^2)/sigma^2));
%                 
% figure; imagesc(gaborEvenFT);
% figure; mesh(gaborEvenFT);
% angle = atan2(X,Y);
% Theta = angle .* (angle >= 0) - (angle .* (angle < 0));
% cosTheta = cos(Theta);
%                 
% % create function to be optimized. IE the l2 distance
% 
% fun = @(x) sum(sum((gaborEvenFT - (exp(-(R/(n/x(1))).^2)) .* ...
%     (polyval([x(3), x(2), 1], R)) .* ...
%     (polyval([x(6), 0, x(5), 0, x(4)], cosTheta))).^2, 2),1)
% 
% xinit = [5, 2, 2, 0.7, 1, 1.3];
% options = optimset('PlotFcns',@optimplotfval);
% [x, fval] = fminsearch(fun, xinit, options)
% current_x = x;
% current_fval = fval;
% 
% 
% x1range = [3:8];
% polyxvalrange = [0.2:0.4:3];
% 
% for idx1 = 1:numel(x1range)
%     x1 = polyxvalrange(idx1);
%     for idx2 = 1:numel(polyxvalrange)
%         x2 = polyxvalrange(idx2);
%         for idx3 = 1:numel(polyxvalrange)
%             x3 = polyxvalrange(idx3);
%             for idx4 = 1:numel(polyxvalrange)
%                 x4 = polyxvalrange(idx4);
%                 for idx5 = 1:numel(polyxvalrange)
%                     x5 = polyxvalrange(idx5);
%                     for idx6 = 1:numel(polyxvalrange)
%                         x6 = polyxvalrange(idx6);
%                         
%                         xinit = [x1, x2, x3, x4, x5, x6];
%                         [x, fval] = fminsearch(fun,xinit,options);
%                         if fval < current_fval
%                             current_fval = fval; 
%                             current_x = x;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

% x = current_x
% approximator =(exp(-R/(n/x(1))).^2) .* polyval([x(3),x(2),1],R) .* ...
%                 (polyval([x(6), 0 , x(5), 0, x(4)], cosTheta));
    

% new search

 n = 100; % size of image
[X, Y] = meshgrid(-n:n);
R = sqrt(X.^2 + Y.^2);

sigma = n/6; 
y0 = n/3;
gaborEvenFT = 0.5 * (exp( -((X).^2+(Y-y0).^2)/sigma^2) + ...
                    exp(-((X).^2+(Y+y0).^2)/sigma^2));
gaborOddFT = 0.5 *   (exp( -((X).^2+(Y-y0).^2)/sigma^2) - ...
                    exp(-((X).^2+(Y+y0).^2)/sigma^2));              
                
figure; imagesc(gaborOddFT);
figure; mesh(gaborOddFT);
angle = atan2(X,Y);
Theta = angle .* (angle >= 0) - (angle .* (angle < 0));
cosTheta = cos(Theta);
                
% create function to be optimized. IE the l2 distance

fun = @(x) sum(sum((gaborOddFT - (exp(-(R/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], R)) .* ...
    (polyval([x(6), 0, x(5), 0], cosTheta))).^2, 2),1)

xinit = [5, 2, 2, 0.7, 1, 1.3];
% options = optimset('PlotFcns',@optimplotfval);
[x, fval] = fminsearch(fun, xinit)
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
approximator =(exp(-R/(n/x(1))).^2) .* polyval([x(4),x(3),x(2)],R) .* ...
                (polyval([x(6), 0 , x(5), 0], cosTheta));

% % start here

n = 100; % size of image
[X, Y] = meshgrid(-n:n);
R = sqrt(X.^2 + Y.^2);

sigma = n/6; 
y0 = n/3;
gaborEvenFT = 0.5 * (exp( -((X).^2+(Y-y0).^2)/sigma^2) + ...
                    exp(-((X).^2+(Y+y0).^2)/sigma^2));
gaborOddFT = 0.5 *   (exp( -((X).^2+(Y-y0).^2)/sigma^2) - ...
                    exp(-((X).^2+(Y+y0).^2)/sigma^2));              
                
figure; imagesc(gaborOddFT);
figure; mesh(gaborOddFT);
angle = atan2(X,Y);
Theta = angle .* (angle >= 0) - (angle .* (angle < 0));
cosTheta = cos(Theta);

fun1 = @(x) sum(sum((gaborOddFT - (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([-0.2668, 0, 0.0852, 0], cosTheta))).^2, 2),1)

xinit = [n/6, 1, 0, 0];
[x,fval] = fminsearch(fun1, xinit);

approximator = (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([-0.2668, 0, 0.0852, 0], cosTheta));

a=x;
fun2 = @(x) sum(sum((gaborOddFT - (exp(-((R-y0)/(n/a(1))).^2)) .* ...
    (polyval([a(4), a(3), a(2)], (R-y0))) .* ...
    (polyval([x(2), 0, x(1), 0], cosTheta))).^2, 2),1)

xinit = [0,0]
[x,fval] = fminsearch(fun2, xinit);
b=x;

approximator = (exp(-((R-y0)/(n/a(1))).^2)) .* ...
    (polyval([a(4), a(3), a(2)], (R-y0))) .* ...
    (polyval([b(2), 0, b(1), 0], cosTheta));

fun1 =@(x) sum(sum((gaborOddFT - (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([b(2), 0, b(1), 0], cosTheta))).^2, 2),1)
xinit = [n/6,1,0,0];

[x,fval] = fminsearch(fun1, xinit);

x = [x,b];
approximator_vals = x;
approximator = (exp(-((R-y0)/(n/x(1))).^2)) .* ...
                polyval([x(4),x(3),x(2)],(R-y0)) .* ...
                (polyval([x(6), 0 , x(5), 0], cosTheta));
funcGabor = @(x) sum(sum(((0.5 *   (exp( -((X).^2+(Y-y0).^2)/x(1)^2) - ...
                    exp(-((X).^2+(Y+y0).^2)/x(1)^2))) - approximator).^2, 1),2);

xinit = [n/6];
sigma_opt = xinit;
[x,fval] = fminsearch(funcGabor, xinit)

oddGabor = (0.5 * (exp( -((X).^2+(Y-y0).^2)/x(1)^2) - ...
                    exp(-((X).^2+(Y+y0).^2)/x(1)^2)));
                
func4 =@(x) sum(sum((oddGabor - (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([x(6), 0, x(5), 0], cosTheta))).^2, 2),1)

xinit = approximator_vals;
[x, fval] = fminsearch(func4, xinit);

approximator = (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([x(6), 0, x(5), 0], cosTheta));

% at this point, x = [5.5856, -2.5088, 0.0222, -0.0001, 0.0702, -0.2451]