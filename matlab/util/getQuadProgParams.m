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