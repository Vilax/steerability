function D = computeDMatrixNorthPoleCap(Theta,L) 
% COMPUTEDMATRIXNORTHPOLECAP   Get D matrix whose eigenvectors give optimally 
% concentrated functions axially symmetric even functions on a cap around 
% the North Pole
%
% D=COMPUTEDMATRIXNORTHPOLECAP(THETA,L) computes for a set max bandwith L and 
% north cap THETA a matrix whose eigenvector encodes the the axially 
% symmetric even polynomial that is maximally concentrated on a cap around 
% the north pole. The basis is taken to be 1, z^2, z^4,..., z^L, where L is
% assumed to be even, and the eigenvector gives the coefficients of the \
% optimally concentrated function written in the monomial basis above. 
% That is the eigenvector is the argmax of the ratio: 
% int_cap f^2 / int_S2 f^2
% If L is odd then the max bandwith will be L-1, and if that is odd an
% error will be thrown. The size of the cap is specified by THETA, which is
% the angle made between the boundary of the cap and the positive z-axis. 

    % validate THETA, L
    assert(Theta>=0, 'k invalid')
    
    L = validateBandwidth(L);

    % compute size of matrix 
    N = L/2+1;
    
    a = [0:N-1]';
    A = 2*repmat(a, [1, N]);
    b = [1:2:L+1];
    B = repmat(b,[N,1]);
    C = A+B; 
    height = cos(Theta);
    cosMat = 1- height .^ C;
    D = cosMat ./ C;
end

function L = validateBandwidth(L)
% perform validation for bandwidth/order L. First ensures it is even, and
% returns error if modified version is less than 0.
    if mod(L,2)~=0
        L=L-1;
    end
    assert(L>=0, 'L invalid');
end