function D = computeNormalizedDMatrixPoleCap(Theta,L) 
% COMPUTENORMALIZEDDMATRIXPOLECAP   Get D matrix whose max eigenvector is 
% the poles-symmetric polynomial optimally concentrated at the pole caps 
%
% D=COMPUTENORMALIZEDDMATRIXPOLECAP(THETA,L) computes for a set max 
% bandwith L and pole cap with sized specified by positive z-axis angle
% THETA a matrix whose eigenvector encodes the the axially 
% symmetric even polynomial that is maximally concentrated on a cap around 
% the (north) pole. The basis is taken to be scaled multiples of the 
% monomials 1, z^2, z^4,..., z^L, where L is
% assumed to be even, (such that the squared norm of a basis vector over a 
% half sphere is 1)and the eigenvector gives the coefficients of the \
% optimally concentrated function written in the monomial basis above. 
% That is the eigenvector is the argmax of the ratio: 
% int_cap f^2 / int_S2 f^2
% If L is odd then the max bandwith will be L-1, and if that is negative an
% error will be thrown. The size of the cap is specified by THETA, which is
% the angle made between the boundary of the cap and the positive z-axis. 

    % validate THETA, L
    assert(Theta>=0, 'k invalid')
    
    L = validateBandwidth(L);

    % compute size of matrix 
    N = L/2+1;
    incr = [1:N]';
    incrRowMat = repmat(incr, [1, N]);
    incrColMat = incrRowMat';
    
    rowScaleMat = sqrt(4 * incrRowMat - 3);
    colScaleMat = rowScaleMat';
    
    % this is the denominator matrix
    denomMat = 2*(incrRowMat + incrColMat) - 3;
    
    % this is coincidentally the matrix for when Theta=pi/2
    B =  rowScaleMat .* colScaleMat ./ denomMat;
    zval = cos(Theta);
    capScaleMat = (1 - zval.^denomMat);
    
    D = B .* capScaleMat;
end

function L = validateBandwidth(L)
% perform validation for bandwidth/order L. First ensures it is even, and
% returns error if modified version is less than 0.
    if mod(L,2)~=0
        L=L-1;
    end
    assert(L>=0, 'L invalid');
end