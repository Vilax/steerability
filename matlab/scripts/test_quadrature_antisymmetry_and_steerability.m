close all; clear all;

% with this script we form special polynomials to steer with, optimizing
% for closeness to quadrature

N_even = 4;
N_odd = 3;
maxDeg = 4; 
numSamples = 400;
[f_even, f_odd, u_even, u_odd, phi,alpha]=makeQuadPoly(maxDeg, numSamples);

f0 = 10;
r0 = 30;
n = 100;
sigmaDenom = n / (0.368*r0);
filtEven = makeSteerFilt3D(n, r0, sigmaDenom, 'interpolation', f_even,phi);
filtOdd = makeSteerFilt3D(n, r0, sigmaDenom, 'interpolation', f_odd,phi);

diff = abs(filtEven) - abs(filtOdd);
quadError = norm(diff(:),2) / norm(filtEven(:),2)
quadError = norm(diff(:),2) / norm(filtOdd(:),2)

WriteMRC(filtEven, 1, 'filtQuadEven.mrc')
WriteMRC(filtOdd, 1, 'filtQuadOdd.mrc')

M_even = (N_even+1)*(N_even+2)/2;
M_odd = (N_odd+1)*(N_odd+2)/2;

DirCosEven =  [ -0.0981    0.0682    0.9928;
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

DirCosOdd = [0.1431    0.3872    0.9108
    0.1751   -0.3962    0.9013
   -0.6130   -0.0307    0.7895
    0.7911    0.0201    0.6113
   -0.4085    0.7575    0.5093
   -0.3662   -0.7877    0.4954
    0.3626    0.8802    0.3061
    0.4121   -0.8638    0.2898
   -0.9278    0.3651    0.0771
   -0.9056   -0.4186    0.0689];

steerAnglesEven = computeSteerBasisAngles3D(N_even, DirCosEven);
PhiEven = inv(steerAnglesEven);
steerAnglesOdd = computeSteerBasisAngles3D(N_odd, DirCosOdd);
PhiOdd = inv(steerAnglesOdd);

steerFiltHyperVolumeEven = makeSteerBasis3D(filtEven, DirCosEven);
nfiltEven = size(steerFiltHyperVolumeEven, 4);

steerFiltHyperVolumeOdd = makeSteerBasis3D(filtOdd, DirCosOdd);
nfiltOdd = size(steerFiltHyperVolumeOdd, 4);

symAxes =[[1, 2,3];
            [3,4,5]];
        
for isym = 1:size(symAxes,1)
    axis = symAxes(isym, :);
    
    dirAnglesEven = getDirectionAngles(N_even, axis);
    dirAnglesOdd = getDirectionAngles(N_odd, axis);
    
    k_even = reshape(PhiEven * dirAnglesEven, [1,1,1,nfiltEven]);
    k_odd = reshape(PhiOdd * dirAnglesOdd, [1,1,1,nfiltOdd]);
    
    steeredFiltEven = sum( k_even .* steerFiltHyperVolumeEven, 4);
    bruteForceFilterEven = rotateFilt3D(filtEven, axis);    
    diffEven = steeredFiltEven - bruteForceFilterEven;
    relativeDiffEven = norm(diffEven(:),2) / norm(bruteForceFilterEven(:),2)
    
    steeredFiltOdd = sum( k_odd .* steerFiltHyperVolumeOdd, 4);
    bruteForceFilterOdd = rotateFilt3D(filtOdd, axis);    
    diffOdd = steeredFiltOdd - bruteForceFilterOdd;
    relativeDiffOdd = norm(diffOdd(:),2) / norm(bruteForceFilterOdd(:),2)
    
    diff = abs(steeredFiltEven) - abs(steeredFiltOdd);
    quadError = norm(diff(:),2) / norm(steeredFiltEven(:),2)
    quadError = norm(diff(:),2) / norm(steeredFiltOdd(:),2)
    
    assert(relativeDiffEven < 0.005)
    assert(relativeDiffOdd < 0.005)
    
end
    