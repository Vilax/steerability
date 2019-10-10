%
%
% testing steerability 2d
%
%
% set up parameters
N = 6;
capSize = pi/4;

Theta = capSize;

% create polynomial
% a = concentratedPoly2D(Theta, N);
% p = even2Poly(a);
% p = p / polyval(p,1);
[f,u, bCos, theta] =   steer2dGeneral(Theta,N);
normvals = computeNormalizationConstant2D(N);
u = u ./ normvals;
p = flipud(u);
p = p/polyval(p,1);
% p = p/ polyval(p,1);

xmin = -200;
xmax = 200;
ymin = -200;
ymax = 200;
bigImgLims = [xmin, xmax, ymin, ymax];
% filt = poly2filter(p, bigImgLims, 0.04,0,5);
r0 = 20;
filt = makeSteerFilt(100, f, theta, r0, 5);
% filt(101:end,:) = -filt(101:end,:);

nonzeroCoefficients = [-N:2:N];
nonzeroCoeffBool = index2boolean(nonzeroCoefficients, N);
M = validateNonzeroCoefficients(nonzeroCoeffBool, N);

steerBasisMat = computeSteerBasis2D(N, nonzeroCoeffBool);
Phi = inv(steerBasisMat);

steerFiltVol = makeSteerBasis(filt, M);
% visualize basis filters:
% nfilt = size(steerFiltVol,3);
% for ifilt = 1:nfilt
%     basisFilt = steerFiltVol(:,:,ifilt);
%     figure; imagesc(basisFilt);
%     figure; mesh(basisFilt);
% end

testThetaArr = [pi/5.5, 3.7*pi/5, 4.7*pi/4, 5*pi/3];
ntestTheta = numel(testThetaArr);
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

% visualizing steering function
thetaRange = [0:0.1:3*pi];
nrotations = numel(thetaRange);
steerCoefficientsArray = zeros(size(thetaRange));
steerCoefficientsArray2 = zeros(size(thetaRange));
for irotation = 1:nrotations
    theta = thetaRange(irotation);
    dirVals = computeSteerDirection2D(N, theta, nonzeroCoeffBool);
    k = real(Phi*dirVals);
    steerCoefficientsArray(irotation) = k(1);
    steerCoefficientsArray2(irotation) = k(2);
end

figure; plot(thetaRange(:), steerCoefficientsArray(:))
figure; plot(thetaRange(:), steerCoefficientsArray2(:))



