function lambda = computeLagrangeCoefficients(H, B, c)
    N = size(H,1);
 
    if nargin < 3
        c = [-1; 1];
    end
    if nargin < 2      
        B = zeros(2,N);
        oddId = [1:2:N];
        evenId = [2:2:N];
        B(1,oddId) = 1;
        B(2,evenId) = 1;
    end
    [he,ho,heo,hoe] = separateMatrixParity(inv(H));
    
    she = sum(he(:));
    sho = sum(ho(:));
    sheo = sum(heo(:));
    shoe = sum(hoe(:));
    
    G = [sho, shoe; sheo, she]; 
%     G = B * inv(H) * B';
    lambda = G \ c;
    lambda1 = lambda(1);
    lambda2 = lambda(2);
end

