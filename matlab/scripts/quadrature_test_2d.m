close all; clear all;

% with this script we form special polynomials to steer with, optimizing
% for closeness to quadrature
capSize = pi/6;
N_even = 8;
N_odd = 7;
maxDeg = 4; 
numSamples = 400;

alpha = computeHilbertPair(N_even, capSize);
[f_even, f_odd, u_even, u_odd, phi] = makeFunctionPair(alpha, numSamples);

f0 = 10;
r0 = 30;
n = 100;
sigmaDenom = n / (0.368*r0);

filtEven = makeSteerFilt(n, f_even, phi, r0, sigmaDenom);
filtOdd = makeSteerFilt(n, f_odd, phi, r0, sigmaDenom);

figure; mesh(filtEven)
hold on; mesh(filtOdd)