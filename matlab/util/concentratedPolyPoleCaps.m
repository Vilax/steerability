function a = concentratedPolyPoleCaps(Theta, L)
% CONCENTRATEDPOLYNORTHPOLECAP returns the coefficients of the even
% polynomial optimally concentrated in the pole caps.
%
% a=CONCENTRATEDPOLYNORTHPOLECAP(Theta,L) computes for a set max bandwidth
% L the coefficients of the polynomial in even power terms that maximizes
% concentration in the north pole cap specified by the z-axis angle Theta.
% The basis of monomials is taken to be z^L, ...,z^2,1 where L is
% assumed to be even (L-1) will be used in place of L otherwise), and so the
% coefficients given in the vector a will correspond to that order 
% 
% see POLYVAL

    % Optimally concentrated polynomial is computed using a Rayleigh
    % quotient, with Cholesky decomposition of the matrix for the inner
    % products of the monomials (as they do not form an orthonormal basis).


    % note that until the end, the computations and the matrix
    % manipulations will be performed with a basis of monomials that have
    % are in increasing order of power, and have been normalized to have
    % square norm 1
    D = 2*computeNormalizedDMatrixPoleCap(Theta,L);
    B = 2*computeNormalizedDMatrixPoleCap(pi/2,L);
     
%     C = chol(B)
%     M = inv(C)*D*inv(C')
%     [V, E] = eig(M);
%     v=V(:,end);
%     C*C'
%     C'*C
%     a = flipud(v); 
    
    C = chol(B);
    M = inv(C')*D*inv(C);
    [V, E] = eig(M)
    [~,iMaxEigenValue] = max(diag(E));
    y = V(:,iMaxEigenValue);
    x = C \ y;

    denormVec = sqrt(4*[1:numel(x)]'-3);
    % a now gives coefficients in un-normalized standard basis z^L,
    % z^{L-2}, ..., z^2, 1 etc.
    a = flipud(x) ./ denormVec;
                            
    
end

