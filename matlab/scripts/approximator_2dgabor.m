n = 100; % size of image
[X, Y] = meshgrid(-n:n);
R = sqrt(X.^2 + Y.^2);

sigma = n/4.5; 
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

func = @(x) sum(sum((gaborOddFT - (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([x(6), 0, x(5), 0], cosTheta)))
xinit = [sigma, 1, 0, 0, 0, 0];.^2, 2),1);
[x,fval] = fminsearch(func, xinit)

b = [x(5), x(6)]

fun1 = @(x) sum(sum((gaborOddFT - (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([b(2), 0, b(1), 0], cosTheta))).^2, 2),1)

xinit = [n/6, 1, 0, 0];
[x,fval] = fminsearch(fun1, xinit);

approximator = (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([b(2), 0, b(1), 0], cosTheta));

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
funcGabor = @(x) sum(sum( ((0.5 * (exp( -((X).^2+(Y-y0).^2)/x(1)^2) - ...
              exp(-((X).^2+(Y+y0).^2)/x(1)^2))) - approximator).^2, 1),2);

xinit = [sigma];

[x, fval] = fminsearch(funcGabor, xinit)
sigma_opt = x;

oddGabor = (0.5 * (exp( -((X).^2+(Y-y0).^2)/x(1)^2) - ...
                    exp(-((X).^2+(Y+y0).^2)/x(1)^2)));

func4 =@(x)sum(sum((oddGabor - (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([x(6), 0, x(5), 0], cosTheta)) ).^2, 1),2)

xinit = approximator_vals;
[x, fval] = fminsearch(func4, xinit)

approximator = (exp(-((R-y0)/(n/x(1))).^2)) .* ...
    (polyval([x(4), x(3), x(2)], (R-y0))) .* ...
    (polyval([x(6), 0, x(5), 0], cosTheta));

% at this point, x = [5.5856, -2.5088, 0.0222, -0.0001, 0.0702, -0.2451]