% first we do it in 2d
clear all; close all;
% set up parameters
N = 6;
capSize = pi/4;
nonzeroCoefficients = [-N:2:N];
nonzeroCoeffBool = index2boolean(nonzeroCoefficients, N);
Theta = capSize;

[f,u, bCos, theta] =   steer2dGeneral(Theta,N);
figure; plot(theta, f);

% first verify steerability of a north-pole oriented filter symmetric over
% the x-axis
n=100;
r0 = 40;

filt = makeSteerFilt(100, f, theta, r0, 5);
figure; imagesc(filt)

steerBasisMat = computeSteerBasis2D(N, nonzeroCoeffBool);
Phi = inv(steerBasisMat);
M = N+1;
steerFiltVol = makeSteerBasis(filt, M);
testThetaArr = [pi/5.5, 3.7*pi/5, 4.7*pi/4, 5*pi/3];
ntestTheta = numel(testThetaArr);
relativeErrorArr = zeros([ntestTheta , 1]);
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

    approxError = (bruteForceRotatedFilt - steeredFilt);
    relativeError = norm(approxError(:),2) / ...
                    norm(bruteForceRotatedFilt(:),2)
    relativeErrorArr(itheta) = relativeError;
end

% now make the negative filt;
multiplierMat = ones(size(filt));
nrows = size(multiplierMat, 1);
center = (nrows+1)/2;
multiplierMat(center+1:end,:) = -1;
negFilt = filt .* multiplierMat; 
figure; imagesc(negFilt);

steerBasisMat = computeSteerBasis2D(N, nonzeroCoeffBool);
Phi = inv(steerBasisMat);
M = N+1;
negSteerFiltVol = makeSteerBasis(negFilt, M);
testThetaArr = [pi/5.5, 3.7*pi/5, 4.7*pi/4, 5*pi/3];
ntestTheta = numel(testThetaArr);
relativeErrorNegArr = zeros([ntestTheta , 1]);

for itheta = 1:ntestTheta
    % extract desired rotation
    thetaTest = testThetaArr(itheta);
    thetaTestDeg = thetaTest * 180/pi;
    % compute manually with imrotate
    bruteForceRotatedFilt=imrotate(negFilt, thetaTestDeg, 'bilinear', 'crop');
    
    % compute e^itheta powers
    dirVals = computeSteerDirection2D(N, thetaTest, nonzeroCoeffBool);
    k = real(Phi * dirVals);
    % compute steered filter 
    steerCoeff = reshape(k, [1,1,numel(k)]);

    steeredFilt = sum((negSteerFiltVol .* steerCoeff), 3);
    % display
    figure; imagesc(bruteForceRotatedFilt);
    figure; imagesc(steeredFilt);

    approxError = (bruteForceRotatedFilt - steeredFilt);
    relativeError = norm(approxError(:),2) / ...
                    norm(bruteForceRotatedFilt(:),2)
    relativeErrorNegArr(itheta) = relativeError;            
end