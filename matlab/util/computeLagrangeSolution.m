function a = computeLagrangeSolution(h, B, c)
    N = size(h,1);
    if nargin < 3
        c = [1; -1];
    end
    if nargin < 2
        B = zeros(2,N);
        oddId = [1:2:N];
        evenId = [2:2:N];
        B(1,oddId) = 1;
        B(2,evenId) = 1;
    end
    G = B * h * B';
    lambda = G \ c;
    lambda1 = lambda(1);
    lambda2 = lambda(2);

    a = h * B' * lambda;
end

