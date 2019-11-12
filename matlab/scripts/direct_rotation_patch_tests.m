% this script runs tests of new rotation tests.
close all; clear all;
N = 5;
Theta = pi/3;
M = 6;
[f,u,bCos,theta]=steer2dNorm(Theta,M);

% retrieve polynomial values and basis and eigenvector (for normalized
% basis)
nonzeroCoefficients = [-N:2:N];
nonzeroCoeffBool = index2boolean(nonzeroCoefficients, N);
M = validateNonzeroCoefficients(nonzeroCoeffBool, N);

nonzeroCoeff = [0,1,0,1,0,1];
[f,u, bCos, theta] =   steer2dGeneral(Theta,N,nonzeroCoeff);
r0 = 50;
sigmaDenom = 5;
filt = makeSteerFilt(100, f, theta, r0, sigmaDenom);

% workflow: test on one or two select angles (old vs new); then test
% overall makesteerbasis function. 

testAngleArr = [12:35:412];
ntests = numel(testAngleArr);
for itest = [1:ntests]
    testAngle = testAngleArr(itest);
    testFilt = imrotate(filt, testAngle, 'bilinear', 'crop');
    filtRotated = makeRotatedFilt(100, f, theta, r0, sigmaDenom, testAngle);
    
    figure; imagesc(testFilt);
    figure; imagesc(filtRotated);
end

nonzeroCoefficients = [-N:2:N];
nonzeroCoeffBool = index2boolean(nonzeroCoefficients, N);
M = validateNonzeroCoefficients(nonzeroCoeffBool, N);
steerFiltVol = makeSteerBasis(100,f,theta,r0,sigmaDenom, M);
testThetaArr = [pi/5.5, 3.7*pi/5, 4.7*pi/4, 5*pi/3];
ntestTheta = numel(testThetaArr);

steerBasisMat = computeSteerBasis2D(N, nonzeroCoeffBool);
Phi = inv(steerBasisMat);

for itheta = 1:ntestTheta
    % extract desired rotation
    thetaTest = testThetaArr(itheta);
    thetaTestDeg = thetaTest * 180/pi;
    % compute manually with imrotate
    bruteForceRotatedFilt=imrotate(filt, thetaTestDeg, 'bilinear', 'crop');
    
    % compute e^itheta powers
    dirVals = computeSteerDirection2D(N, thetaTest, nonzeroCoeffBool);
    k = real(Phi * dirVals)
    % compute steered filter 
    steerCoeff = reshape(k, [1,1,numel(k)]);
    steeredFilt = sum((steerFiltVol .* steerCoeff), 3);
    % display
    figure; imagesc(bruteForceRotatedFilt);
    figure; imagesc(steeredFilt);
%     figure; mesh(bruteForceRotatedFilt);
%     figure; mesh(steeredFilt);
    approxError = (bruteForceRotatedFilt - steeredFilt);
    tmpDenom = bruteForceRotatedFilt;
    tmpDenom(tmpDenom==0)=1;
    figure; mesh(approxError)
    
end

%% now for the 3d part

close all; clear all;

% steering 3d filters
% makes polynomials using interpolation method and direct computation

% now for odd filters
close all;
N = 3;
poleCap = pi/3;

nonzeroBool = [0 1 0 1];

[f,u,bCos,phi] = steer3dGeneral(poleCap, N, nonzeroBool);

figure; plot(phi, f);
% phiVals = [0:0.1:pi/2];
% cosVals = cos(phiVals);
% figure; plot(phiVals, polyval(p,cosVals));

f0 = 10;
r0 = 30;
n = 100;
sigmaDenom = 100/(0.368*r0);
filt = makeSteerFilt3D(n, r0, sigmaDenom, 'interpolation', f,phi);

% first test a single rotation

testAxis = [1.4, 2, 0.9];
testFilt = rotateFilt3D(filt, testAxis); % expected
orientedFilt = makeOrientedFilter(f, phi, testAxis, 'gaussian', ...
                                                        n,r0,sigmaDenom);

diff = testFilt - orientedFilt;
relDiff = norm(diff(:),2) / norm(testFilt(:),2)
sliceTest = squeeze(testFilt(:,101,:));
orientTest = squeeze(orientedFilt(:,101,:));
figure;imagesc(sliceTest);
figure; imagesc(orientTest);

%% testing it for a steering task

close all;
N = 3;
poleCap = pi/3;

nonzeroBool = [0 1 0 1];

[f,u,bCos,phi] = steer3dGeneral(poleCap, N, nonzeroBool);

figure; plot(phi, f);
% phiVals = [0:0.1:pi/2];
% cosVals = cos(phiVals);
% figure; plot(phiVals, polyval(p,cosVals));

f0 = 10;
sigmaDenom = 5;
r0 = 30;
n = 100;
filt = makeSteerFilt3D(n, r0, sigmaDenom, 'interpolation', f,phi);

M = (N+1)*(N+2)/2;
% [V,~,~,~] = ParticleSampleSphere('N', 2*M);
% DirCos = V(1:M, :);
DirCos = [0.1431    0.3872    0.9108
    0.1751   -0.3962    0.9013
   -0.6130   -0.0307    0.7895
    0.7911    0.0201    0.6113
   -0.4085    0.7575    0.5093
   -0.3662   -0.7877    0.4954
    0.3626    0.8802    0.3061
    0.4121   -0.8638    0.2898
   -0.9278    0.3651    0.0771
   -0.9056   -0.4186    0.0689];

steerAngles = computeSteerBasisAngles3D(N, DirCos);
[Q,R] = qr(steerAngles);
Phi = inv(steerAngles);

steerFiltHyperVolume = makeSteerBasis3D(f,phi,DirCos,'gaussian',n,r0,sigmaDenom);
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
    
    diff = steeredFilt - bruteForceFilter;
    normalizedDiff = norm(diff(:),2) / norm(steeredFilt(:),2)
    assert(normalizedDiff < 0.005);
end
