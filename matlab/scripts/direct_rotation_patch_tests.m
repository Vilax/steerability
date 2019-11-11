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
