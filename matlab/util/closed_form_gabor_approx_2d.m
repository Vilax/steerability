% first run gabor approximations
n=100;
[x,approximator,gaborSigma,gaborFT]=approxGaborBCD(n,6,4,10,'even');

assert(numel(x) == 6)
b = [x(5), x(6)];
angularPolyCoeff = [b(2), 0, b(1), 0];

[X,Y] = meshgrid(-n:n);
R = sqrt(X.^2+Y.^2);
angle = atan2(X,Y);
Theta = angle .* (angle >= 0) - (angle .* (angle < 0));
cosTheta = cos(Theta);
angularPart = (polyval([x(6), 0, x(5), 0, 0], cosTheta));

rmax = floor(max(R(:)));
degMax = 360;
imP = ImToPolar(mat2gray(gaborFT), 0, 1, rmax+1, degMax);

angleDeg = [0:359];
theta = angleDeg * pi / 180;
costhetaq = cos(theta);

pvals = polyval(angularPolyCoeff, costhetaq);

% now multiply pvals with imP in the sense of pvals(id) times the column of
% imP, then sum imP to collapse into a single col. This gives h(r) from 0
% to rmax. Then extend h(r) to a 2d map if desired.

colProduct = imP * (diag(pvals)); 
intProduct = sum(imP, 2);

radialPartPolar = intProduct / (norm(angularPart(:),2)^2);

Rq = reshape(R, [numel(R(:)), 1]); 
radialCart = interp1([0:rmax]', radialPartPolar, Rq); 

radialPart = reshape(radialCart, size(R));

approximator = radialPart .* angularPart;

figure; plot([0:rmax]', radialPartPolar)
figure; imagesc(angularPart)
figure; imagesc(radialPart)
figure; imagesc(approximator)
figure; contour(approximator,40)
figure; mesh(approximator)

% now compute the bias of the maximum

gaborMax = max(gaborFT(:));
[j1,j2] = ind2sub(size(gaborFT), find(gaborFT == gaborMax));
ygaborMax = j1

approxMax = max(approximator(:));
[i1,i2] = ind2sub(size(approximator), find(approximator==approxMax));
yapproxMax = i1

bias = abs(j1-i1)
