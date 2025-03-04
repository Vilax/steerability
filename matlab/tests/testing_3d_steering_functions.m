% test rotation 3d
close all;
N = 4;
northPole = [0,0,1];
poleCap = pi/4;

a = concentratedPolyPoleCaps(poleCap, N);
p = even2Poly(a);
p = p/polyval(p,1);


phiVals = [0:0.01:pi/2];
cosVals = cos(phiVals);
figure; plot(phiVals, polyval(p, cosVals));

f0 = 10;
sigmaDenom = 5;
r0=40;
n=100;
filt = makeSteerFilt3D(n, r0, sigmaDenom, 'poly', p);

axis1 = [3 1 2];
filt1bf = rotateFilt3D(filt, axis1);
filtneutralbf = rotateFilt3D(filt1bf, [0,0,1], [3 1 2]);

diff = filt - filtneutralbf;
normalizedDiff = norm(diff(:),2) / norm(filt(:),2)
assert(normalizedDiff < 0.005);
% voxelDiffNormalized = diff ./ filt;
% normDiffVoxelwise = norm(voxelDiffNormalized(:), 2)
% assert(normDiffVoxelwise < 0.005);



% testing that steer3d general and explicit dmatrix methods in 3d give
% similar polynomials

close all;
N = 4;
northPole = [0,0,1];
poleCap = pi/4;

a = concentratedPolyPoleCaps(poleCap, N);
p = even2Poly(a);
p = p/polyval(p,1);

[f,u,bCos,phi] = steer3dGeneral(poleCap, N);
cosVals = cos(phi);
figure; plot(phi, f); hold on;
plot(phi, polyval(p, cosVals));


