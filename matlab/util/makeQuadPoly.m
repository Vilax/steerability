function [f1, f2, u1, u2, phi, alpha] = makeQuadPoly(N, numSamples)

    DEFAULT_NUMSAMPLES = 300;
    if nargin < 2
        numSamples = DEFAULT_NUMSAMPLES;
    end
    
    [H, f, A, b, Aeq, beq] = getQuadProgParams(N);
    alpha = quadprog(H,f,A,b,Aeq,beq);
    mismatch=sqrt(alpha'*H*alpha)
    
    u1 = zeros(size(alpha));
    u2 = zeros(size(alpha));
    u1(1:2:N) = alpha(1:2:N);
    u2(2:2:N) = alpha(2:2:N);
    u2 = [u2(2:end); 0] / sum(u2(:));
    u1 = [u1;0] / sum(u2(:));
    
    dphi = pi/(numSamples - 1);
    phi = [0:dphi:pi]';
    cosPhi = cos(phi);
    
    f1 = polyval(u1, cosPhi);
    f2 = polyval(u2, cosPhi);
    %     figure; plot(cosPhi, f1);
%     hold on; plot(cosPhi,f2);

    
    x = [0:1/(numSamples-1):1]';
    g1 = polyval(u1, x);
    g2 = polyval(u2, x);
%     figure; plot(x, f1); hold on;
%     plot(x, f2);
end

function [H, f, A, b, Aeq, beq] = getQuadProgParams(N)
    H = getQuadProgObjMatrix(N);
    f = zeros(N,1);
    A = [];
    b = [];
    Aeq=zeros(2,N);
    for i=1:2:N
        Aeq(1,i)=1;
    end
    for i=2:2:N
        Aeq(2,i)=1;
    end
    beq = [1,-1];
    if mod(N,2) ~= 0
        beq = -1 * beq;
    end
end
