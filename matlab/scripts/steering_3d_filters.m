% steering 3d filters
% makes polynomials using interpolation method and direct computation
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

M = (N+1)*(N+2)/2;
% [V,~,~,~] = ParticleSampleSphere('N', 2*M);
% DirCos = V(1:M, :);
DirCos =  [ -0.0981    0.0682    0.9928;
    0.4509   -0.2945    0.8426;
   -0.2773   -0.5078    0.8156;
    0.4556    0.3731    0.8082;
   -0.6677    0.2072    0.7150;
   -0.1407    0.6864    0.7135;
    0.2475   -0.8244    0.5090;
   -0.8258   -0.3546    0.4385;
    0.8952    0.1935    0.4016;
    0.8089   -0.4851    0.3321;
    0.4799    0.8158    0.3228;
   -0.4383   -0.8555    0.2756;
   -0.6724    0.6930    0.2602;
   -0.1040    0.9883    0.1115;
   -0.9801    0.1936    0.0435];

steerAngles = computeSteerBasisAngles3D(N, DirCos);
[Q,R] = qr(steerAngles);
Phi = inv(steerAngles);

steerFiltHyperVolume = makeSteerBasis3D(filt, DirCos);
nfilt = size(steerFiltHyperVolume, 4);

symAxes =[[1, 2,3];
            [3,4,5]];
for isym = 1:size(symAxes,1)        
    axis = symAxes(isym, :);
    alpha = axis(1) /norm(axis);
    beta = axis(2)/norm(axis);
    gamma = axis(3)/norm(axis);
    
    powers = dirCosPowers3D(N);
    alphaPowers = powers(:,1);
    betaPowers = powers(:,2);
    gammaPowers = powers(:,3);
    
    dirAlpha = (alpha .^ alphaPowers);
    dirBeta = (beta .^ betaPowers);
    dirGamma = (gamma .^ gammaPowers);
    dirAngles = dirAlpha .* dirBeta .* dirGamma;
    

    k = steerAngles \ dirAngles;
%     k = R \ (Q'*dirAngles);
    k = reshape(k,[1,1,1, nfilt]);    
    
    steeredFilt = sum( k .* steerFiltHyperVolume, 4);
    bruteForceFilter = rotateFilt3D(filt, axis);
end

% now for odd filters
% close all;
% N = 3;
% poleCap = pi/3;
% 
% nonzeroBool = [0 1 0 1];
% 
% % HERE DEFINE AND NORMALIZE ODD POLYNOMIAL P
% 
% phiVals = [0:0.1:pi/2];
% cosVals = cos(phiVals);
% figure; plot(phiVals, polyval(p,cosVals));
% 
% f0 = 10;
% sigmaDenom = 5;
% r0 = 30;
% n = 100;
% filt = makeSteerFilt3D(n, r0, sigmaDenom, 'poly', p);


